//! 1:1 translation of OCCT `ShapeFix_Face` — part C (the impl submodule of
//! [`super::face_a::ShapeFixFace`], the OCCT continued-file convention):
//! `FixMissingSeam` (cxx L1722-2326), `FixSmallAreaWire` (cxx L2331-2394),
//! `FixLoopWire` (cxx L2478-2642), `SplitEdge` x2 (cxx L2646-2817),
//! `FixIntersectingWires` (cxx L2821-2825), `FixWiresTwoCoincEdges` (cxx
//! L2829-2901), `FixSplitFace` (cxx L2905-3010), `FixPeriodicDegenerated`
//! (cxx L3101-3259) and the file statics `CheckWire` (cxx L1652-1718),
//! `FindNext` (cxx L2398-2456), `isClosed2D` (cxx L2458-2474),
//! `IsPeriodicConicalLoop` (cxx L3018-3097).
//!
//! GAP carriers (the iron rule: dependency + anchor + OCCT failure path):
//! - [`ShapeFixComposeShellGap`] — OCCT `ShapeFix_ComposeShell`
//!   (ShapeFix_ComposeShell.cxx, 3,606 LOC; docket C-on-demand row, zero
//!   oracle in the 36-case face) reduced to the
//!   Init/ClosedMode/SetContext/SetMaxTolerance/Perform/Result call shape of
//!   `FixMissingSeam` (cxx L2252-2266); `Perform` keeps OCCT's "nothing
//!   composed" path (the result is the input face).
//! - `BRepBuilderAPI_MakeFace(surf, tol)` — the natural-bound wire
//!   construction (BRepLib_MakeFace) is kernel scope; the face is created
//!   with the natural-restriction flag and no wires (annotated at the
//!   [`super::face_b`] sites).

use std::collections::HashMap;
use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Line2d, Surface3};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType, TShape};

use crate::shhealing::shape_analysis::analysis;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_analysis::wire::ShapeAnalysisWire;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes, topexp_explorer};
use crate::shhealing::shape_build::edge::{builder_range_on_face, ShapeBuildEdge};
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::composite_surface::ShapeExtendCompositeSurface;
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::face_a::{occt_is_same, WireListMap};
use crate::shhealing::shape_fix::face_b::set_edge_degenerated;
use crate::shhealing::shape_fix::intersection_tool::ShapeFixIntersectionTool;
use crate::shhealing::shape_fix::split_tool::{
    brep_tool_surface_loc, topexp_vertices, ShapeFixSplitTool,
};
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::shhealing::shape_fix::wire::{
    brep_tool_degenerated, shape_oriented, ShapeFixWire,
};
use crate::topalgo::brep_top_adaptor::fclass2d_topol::FClass2dTopol;

use super::face_a::ShapeFixFace;
use super::intersection_tool::{bnd_lib_add2d_curve, shape_is_equal};

impl ShapeFixFace {
    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L1722-2326 — FixMissingSeam.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixMissingSeam (cxx L1722-2326): detects and fixes
    /// the special case when a face on a closed surface is given by two
    /// wires closed in 3d but with a gap in 2d — creates a new wire from the
    /// two and adds the missing seam edge (the :i7 abv 18 Sep 98 algorithm).
    pub fn fix_missing_seam(&mut self, brep: &mut BRep) -> bool {
        // OCCT L1724-1727.
        if self.my_surf.is_none() {
            return false;
        }

        // OCCT L1729-1735.
        let (uclosed, vclosed) = {
            let sas = self.my_surf.as_mut().unwrap();
            (sas.is_u_closed(CONFUSION), sas.is_v_closed(CONFUSION))
        };
        if !uclosed && !vclosed {
            return false;
        }

        // OCCT L1737-1741.
        if self.base.my_context.is_some() {
            let f = self.my_face.clone();
            let s = self.context_apply(brep, &f);
            self.my_face = s;
        }

        // OCCT L1743-1751: %pdn — the surface should be made periodic before
        // (see ShapeCustom_Surface)!
        {
            let sas = self.my_surf.as_ref().unwrap();
            if let Surface3::BSpline(bspl) = sas.surface() {
                let bspl = bspl.clone();
                let u_per = rcad_kernel::geom::SurfaceEval::is_u_periodic(&Surface3::BSpline(bspl.clone()));
                let v_per = rcad_kernel::geom::SurfaceEval::is_v_periodic(&Surface3::BSpline(bspl));
                if !u_per && !v_per {
                    return false;
                }
            }
        }

        // OCCT L1753-1756.
        let (mut su_f, mut su_l, mut sv_f, mut sv_l) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        self.my_surf.as_ref().unwrap().bounds(&mut su_f, &mut su_l, &mut sv_f, &mut sv_l);
        let (mut f_u1, mut f_u2, mut f_v1, mut f_v2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        {
            // OCCT L1756: BRepTools::UVBounds(myFace, fU1, fU2, fV1, fV2) —
            // the pcurve bounds (the GetFaceUVBounds service).
            analysis::get_face_uv_bounds(brep, &self.my_face, &mut f_u1, &mut f_u2, &mut f_v1, &mut f_v2);
        }

        // OCCT L1758-1802: pdn OCC55 — faces without the wires.
        if su_f.abs() >= rcad_kernel::precision::INFINITE_VALUE
            || su_l.abs() >= rcad_kernel::precision::INFINITE_VALUE
        {
            if su_f.abs() >= rcad_kernel::precision::INFINITE_VALUE {
                su_f = f_u1;
            }
            if su_l.abs() >= rcad_kernel::precision::INFINITE_VALUE {
                su_l = f_u2;
            }
            if (su_l - su_f).abs() < PCONFUSION {
                if su_f.abs() >= rcad_kernel::precision::INFINITE_VALUE {
                    su_f -= 1000.0;
                } else {
                    su_l += 1000.0;
                }
            }
        }
        if sv_f.abs() >= rcad_kernel::precision::INFINITE_VALUE
            || sv_l.abs() >= rcad_kernel::precision::INFINITE_VALUE
        {
            if sv_f.abs() >= rcad_kernel::precision::INFINITE_VALUE {
                sv_f = f_v1;
            }
            if sv_l.abs() >= rcad_kernel::precision::INFINITE_VALUE {
                sv_l = f_v2;
            }
            if (sv_l - sv_f).abs() < PCONFUSION {
                if sv_f.abs() >= rcad_kernel::precision::INFINITE_VALUE {
                    sv_f -= 1000.0;
                } else {
                    sv_l += 1000.0;
                }
            }
        }

        // OCCT L1804-1808.
        let u_range = (su_l - su_f).abs().min(rcad_kernel::precision::INFINITE_VALUE);
        let v_range = (sv_l - sv_f).abs().min(rcad_kernel::precision::INFINITE_VALUE);
        let mut ismodeu = 0i32;
        let mut ismodev = 0i32; // szv#4:S4163:12Mar99 was Boolean
        let mut isdeg1 = 0i32;
        let mut isdeg2 = 0i32;

        // OCCT L1810-1822.
        let mut ws: Vec<Shape> = Vec::new();
        let mut a_seq_non_manif: Vec<Shape> = Vec::new();
        for value in iter_subshapes(brep, &self.my_face, false, true) {
            if value.shape_type() != ShapeType::Wire
                || (value.orientation != Orientation::Forward
                    && value.orientation != Orientation::Reversed)
            {
                a_seq_non_manif.push(value);
                continue;
            }
            ws.push(value);
        }

        // OCCT L1824-1881.
        let mut w1 = Shape::null();
        let mut w2 = Shape::null();
        let mut i = 1i32;
        while i <= ws.len() as i32 {
            let wire = ws[(i - 1) as usize].clone();
            let (mut isuopen, mut isvopen) = (0i32, 0i32);
            let mut isdeg = false;
            if !check_wire(brep, &wire, &self.my_face, u_range, v_range, &mut isuopen, &mut isvopen, &mut isdeg) {
                i += 1;
                continue;
            }
            if w1.is_null() {
                w1 = wire;
                ismodeu = isuopen;
                ismodev = isvopen;
                isdeg1 = if isdeg { i } else { 0 };
            } else if w2.is_null() {
                if ismodeu == -isuopen && ismodev == -isvopen {
                    w2 = wire;
                    isdeg2 = if isdeg { i } else { 0 };
                } else if ismodeu == isuopen && ismodev == isvopen {
                    w2 = wire.clone();
                    isdeg2 = if isdeg { i } else { 0 };
                    // :abv 29.08.01: if wires are contraversal, reverse one of
                    // them; if the first one is a single degenerated edge,
                    // reverse it; else the second
                    if isdeg1 != 0 {
                        w1.orientation = match w1.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            o => o,
                        };
                        ismodeu = -ismodeu;
                        ismodev = -ismodev;
                    } else {
                        w2.orientation = match w2.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            o => o,
                        };
                    }
                }
            }
            //    else return false; //  abort
            else {
                // :abv 30.08.09: if more than one open wires and more than
                // two of them are completely degenerated, remove any of them
                if isdeg || isdeg1 != 0 || isdeg2 != 0 {
                    let remove_idx = if isdeg { i } else if isdeg2 != 0 { isdeg2 } else { isdeg1 };
                    ws.remove((remove_idx - 1) as usize);
                    w1 = Shape::null();
                    w2 = Shape::null();
                    // OCCT L1877: i = 0; continue — the scan restarts.
                    i = 0;
                }
            }
            i += 1;
        }

        // OCCT L1883-1893.
        let a_tor_surf = match self.my_surf.as_ref().unwrap().surface() {
            Surface3::Torus(t) => Some(t.clone()),
            _ => None,
        };
        let mut an_is_degenerated_tor = match &a_tor_surf {
            Some(t) => t.major_radius < t.minor_radius,
            None => false,
        };
        // if the second wire is not null, we don't need mark the torus as
        // degenerated and should process it as a regular one.
        if an_is_degenerated_tor && !w2.is_null() {
            an_is_degenerated_tor = false;
        }

        // OCCT L1895-1992.
        if w1.is_null() {
            return false;
        } else if w2.is_null() {
            // For spheres and BSpline cone-like surfaces (bug 24055).
            let mut p = DVec2::ZERO;
            let mut d = DVec2::ZERO;
            let a_range;

            if ismodeu != 0 && an_is_degenerated_tor {
                // OCCT L1909-1919.
                let a_ra = a_tor_surf.as_ref().unwrap().major_radius;
                let a_ri = a_tor_surf.as_ref().unwrap().minor_radius;
                let a_phi = (-a_ra / a_ri).acos();
                p = DVec2::new(0.0, if ismodeu > 0 { std::f64::consts::PI + a_phi } else { a_phi });

                let a_x_coord = -(ismodeu as f64);
                d = DVec2::new(a_x_coord, 0.0);
                a_range = 2.0 * std::f64::consts::PI;
            } else if ismodeu != 0
                && matches!(self.my_surf.as_ref().unwrap().surface(), Surface3::Sphere(_))
            {
                // OCCT L1920-1926.
                p = DVec2::new(
                    if ismodeu < 0 { 0.0 } else { 2.0 * std::f64::consts::PI },
                    ismodeu as f64 * 0.5 * std::f64::consts::PI,
                );
                let a_x_coord = -(ismodeu as f64);
                d = DVec2::new(a_x_coord, 0.0);
                a_range = 2.0 * std::f64::consts::PI;
            } else if ismodev != 0
                && matches!(self.my_surf.as_ref().unwrap().surface(), Surface3::BSpline(_))
            {
                // OCCT L1927-1948.
                let sas = self.my_surf.as_ref().unwrap();
                let mut u_coord;
                if sas.value(DVec2::new(su_f, sv_f)).distance(sas.value(DVec2::new(su_f, (sv_f + sv_l) / 2.0)))
                    < CONFUSION
                {
                    u_coord = su_f;
                } else if sas.value(DVec2::new(su_l, sv_f)).distance(sas.value(DVec2::new(su_l, (sv_f + sv_l) / 2.0)))
                    < CONFUSION
                {
                    u_coord = su_l;
                } else {
                    return false;
                }
                let _ = &mut u_coord;

                p = DVec2::new(u_coord, if ismodev < 0 { 0.0 } else { v_range });
                d = DVec2::new(0.0, -(ismodev as f64));
                a_range = v_range;
            } else if ismodeu != 0
                && matches!(self.my_surf.as_ref().unwrap().surface(), Surface3::BSpline(_))
            {
                // OCCT L1949-1971.
                let sas = self.my_surf.as_ref().unwrap();
                let mut v_coord;
                if sas.value(DVec2::new(su_f, sv_f)).distance(sas.value(DVec2::new((su_f + su_l) / 2.0, sv_f)))
                    < CONFUSION
                {
                    v_coord = sv_f;
                } else if sas.value(DVec2::new(su_l, sv_l)).distance(sas.value(DVec2::new((su_f + su_l) / 2.0, sv_l)))
                    < CONFUSION
                {
                    v_coord = sv_l;
                } else {
                    return false;
                }
                let _ = &mut v_coord;

                p = DVec2::new(if ismodeu < 0 { 0.0 } else { u_range }, v_coord);
                let a_x_coord = -(ismodeu as f64);
                d = DVec2::new(a_x_coord, 0.0);
                a_range = u_range;
            } else {
                return false;
            }

            // OCCT L1977-1991.
            let line2d = Curve2d::Line(Line2d {
                origin: p,
                direction: d.normalize_or_zero(),
            });
            let mut b = BRepBuilder::new();
            let pnt = self.my_surf.as_ref().unwrap().value(p);
            let v = b.add_vertex(brep, pnt, CONFUSION);
            let v_fwd = shape_oriented(&v, Orientation::Forward);
            let v_rev = shape_oriented(&v, Orientation::Reversed);
            // OCCT L1979-1982: B.MakeEdge(edge); B.Degenerated(edge, true);
            // B.UpdateEdge(edge, line, myFace, Confusion); B.Range(edge,
            // myFace, 0., aRange) — the rcad edge constructor carries the
            // vertices (the B.Add(edge, V) pairs below).
            let mut edge = b.add_edge(brep, None, v_fwd.clone(), v_rev, [0.0, 0.0]);
            set_edge_degenerated(brep, &edge, true);
            b.update_edge_pcurve(brep, edge.clone(), line2d.clone(), self.my_face.clone(), CONFUSION);
            builder_range_on_face(brep, &edge, &self.my_face, 0.0, a_range);
            // OCCT L1983-1988: B.MakeVertex; V.Orientation x2; B.Add(edge, V)
            // — carried by the add_tedge vertex pair.
            w2 = brep.add_twire(Vec::new());
            builder_add(brep, &w2, &edge);
            ws.push(w2.clone());
            let _ = &mut edge;
        }

        // OCCT L1994-2006.
        let mut uf = su_f;
        let mut vf = sv_f;
        let coord = if ismodeu != 0 { 1usize } else { 0usize };
        let isneg = if ismodeu != 0 { ismodeu } else { -ismodev };
        let period = if ismodeu != 0 { u_range } else { v_range };
        let mut m1 = [[0.0f64; 2]; 2];
        let mut m2 = [[0.0f64; 2]; 2];
        {
            let (mut a0, mut b0, mut c0, mut d0) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            let s_empty = brep.empty_copied(&self.my_face);
            let mut s = s_empty;
            builder_add(brep, &s, &w1);
            analysis::get_face_uv_bounds(brep, &s, &mut a0, &mut b0, &mut c0, &mut d0);
            m1[0][0] = a0;
            m1[0][1] = b0;
            m1[1][0] = c0;
            m1[1][1] = d0;
            let (mut a0, mut b0, mut c0, mut d0) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            let s_empty = brep.empty_copied(&self.my_face);
            let mut s = s_empty;
            builder_add(brep, &s, &w2);
            analysis::get_face_uv_bounds(brep, &s, &mut a0, &mut b0, &mut c0, &mut d0);
            m2[0][0] = a0;
            m2[0][1] = b0;
            m2[1][0] = c0;
            m2[1][1] = d0;
        }

        // OCCT L2008-2021.
        if !vclosed || !uclosed || an_is_degenerated_tor {
            let delta_other = 0.5 * (m2[coord][0] + m2[coord][1]) - 0.5 * (m1[coord][0] + m1[coord][1]);
            if delta_other * (isneg as f64) < 0.0 {
                w1.orientation = match w1.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };
                w2.orientation = match w2.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };
            }
        }

        // OCCT L2023-2034: sort original wires.
        let mut sfw = ShapeFixWire::new();
        sfw.set_face_with_surface(brep, &self.my_face, super::face_a::sas_dup(self.my_surf.as_ref().unwrap()));
        sfw.set_precision(self.base.my_precision);
        let wd1 = WireData::new_from_wire(brep, &w1, true, false);
        let wd2 = WireData::new_from_wire(brep, &w2, true, false);
        sfw.load_wire_data(brep, wd1);
        sfw.fix_reorder(brep, false);
        let wd1 = sfw.my_analyzer.my_wire.take().unwrap();
        sfw.load_wire_data(brep, wd2);
        sfw.fix_reorder(brep, false);
        let wd2 = sfw.my_analyzer.my_wire.take().unwrap();
        let w11 = wd1.wire(brep);
        let w21 = wd2.wire(brep);

        // OCCT L2036-2064: :abv 29.08.01 — reconstruct face taking into
        // account reversing.
        let dummy = brep.empty_copied(&self.my_face);
        let mut tmp_f = dummy;
        tmp_f.orientation = Orientation::Forward;
        for i in 1..=(ws.len() as i32) {
            let mut wire = ws[(i - 1) as usize].clone();
            if wire.is_same(&w1) {
                wire = w11.clone();
            } else if wire.is_same(&w2) {
                wire = w21.clone();
            } else {
                // other wires (not boundary) are considered as holes; make
                // sure to have them oriented accordingly
                let curface_empty = brep.empty_copied(&tmp_f);
                let mut curface = curface_empty;
                builder_add(brep, &curface, &wire);
                curface.orientation = self.my_face.orientation;
                if analysis::is_outer_bound(brep, &curface) {
                    wire.orientation = match wire.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        o => o,
                    };
                }
            }
            builder_add(brep, &tmp_f, &wire);
        }

        tmp_f.orientation = self.my_face.orientation;

        // OCCT L2066-2136: a special kind of FixShifted for torus-like
        // surfaces (tr9_r0501-ug.stp #187640).
        if uclosed && vclosed && !an_is_degenerated_tor {
            let shiftw2 = analysis::adjust_by_period(
                0.5 * (m2[coord][0] + m2[coord][1]),
                0.5 * (m1[coord][0] + m1[coord][1] + (isneg as f64) * (period + PCONFUSION)),
                period,
            );
            m1[coord][0] = m1[coord][0].min(m2[coord][0] + shiftw2);
            m1[coord][1] = m1[coord][1].max(m2[coord][1] + shiftw2);
            for value in iter_subshapes(brep, &tmp_f, false, true) {
                if value.shape_type() != ShapeType::Wire {
                    continue;
                }
                let w = value;
                if occt_is_same(&w, &w11) {
                    continue;
                }
                let shift;
                if occt_is_same(&w, &w21) {
                    shift = shiftw2;
                } else {
                    let (mut a0, mut b0, mut c0, mut d0) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                    let s_empty = brep.empty_copied(&tmp_f);
                    let mut s = s_empty;
                    builder_add(brep, &s, &w);
                    analysis::get_face_uv_bounds(brep, &s, &mut a0, &mut b0, &mut c0, &mut d0);
                    m2[0][0] = a0;
                    m2[0][1] = b0;
                    m2[1][0] = c0;
                    m2[1][1] = d0;
                    shift = analysis::adjust_by_period(
                        0.5 * (m2[coord][0] + m2[coord][1]),
                        0.5 * (m1[coord][0] + m1[coord][1]),
                        period,
                    );
                }
                if shift != 0.0 {
                    // OCCT L2106-2119: gp_Vec2d V; V.SetCoord(coord + 1, shift).
                    let sae = ShapeAnalysisEdge::new();
                    for iw in iter_subshapes(brep, &w, true, true) {
                        let e = iw;
                        let mut c: Option<Curve2d> = None;
                        let mut a = 0.0f64;
                        let mut bb = 0.0f64;
                        if !sae.pcurve_face(brep, &e, &tmp_f, &mut c, &mut a, &mut bb, false) {
                            continue;
                        }
                        // OCCT L2118: C->Translate(V).
                        let c = c.unwrap();
                        let offset = if coord == 0 {
                            DVec2::new(shift, 0.0)
                        } else {
                            DVec2::new(0.0, shift)
                        };
                        let c = rcad_kernel::geom::translate_curve2d(&c, offset);
                        let sbe = ShapeBuildEdge;
                        sbe.replace_pcurve(brep, &e, &c, &tmp_f);
                    }
                }
            }
            // OCCT L2122-2135: abv 05 Feb 02: OCC34 — select the proper split
            // place by V to avoid extra intersections.
            if m1[coord][1] - m1[coord][0] <= period {
                let other = 0.5 * (m1[coord][0] + m1[coord][1] - period);
                if ismodeu != 0 {
                    vf = other;
                } else {
                    uf = other;
                }
            }
        }

        // OCCT L2138-2224: find the best place by u and v to insert a seam.
        let sae = ShapeAnalysisEdge::new();
        let mut found_u = 0i32;
        let mut found_v = 0i32;
        let nb1 = wd1.nb_edges();
        let nb2 = wd2.nb_edges();
        let mut i1 = 1i32;
        'seam: while i1 <= nb1 + nb2 {
            let edge1 = if i1 <= nb1 { wd1.edge(i1) } else { wd2.edge(i1 - nb1) };
            let mut c2d: Option<Curve2d> = None;
            let mut f = 0.0f64;
            let mut l = 0.0f64;
            if !sae.pcurve_face(brep, &edge1, &tmp_f, &mut c2d, &mut f, &mut l, true) {
                return false;
            }
            let c2d = c2d.unwrap();
            let mut pos1 = Curve2dEval::point_at(&c2d, l);
            // the best place is end of edge which is nearest to 0
            let mut skip_u = !uclosed;
            if uclosed && ismodeu != 0 {
                pos1.x += analysis::adjust_by_period(pos1.x, su_f, u_range);
                if found_u == 2 && pos1.x.abs() > uf.abs() {
                    skip_u = true;
                } else if found_u == 0 || (found_u == 1 && pos1.x.abs() < uf.abs()) {
                    found_u = 1;
                    uf = pos1.x;
                }
            }
            let mut skip_v = !vclosed;
            if vclosed && ismodeu == 0 {
                pos1.y += analysis::adjust_by_period(pos1.y, sv_f, v_range);
                if found_v == 2 && pos1.y.abs() > vf.abs() {
                    skip_v = true;
                } else if found_v == 0 || (found_v == 1 && pos1.y.abs() < vf.abs()) {
                    found_v = 1;
                    vf = pos1.y;
                }
            }
            if skip_u && skip_v {
                if i1 <= nb1 {
                    i1 += 1;
                    continue 'seam;
                } else {
                    break 'seam;
                }
            }
            // or yet better - if it is end of some edges on both wires
            let mut i2 = 1i32;
            while i1 <= nb1 && i2 <= nb2 {
                let edge2 = wd2.edge(i2);
                let mut c2d2: Option<Curve2d> = None;
                let mut f2 = 0.0f64;
                let mut l2 = 0.0f64;
                if !sae.pcurve_face(brep, &edge2, &tmp_f, &mut c2d2, &mut f2, &mut l2, true) {
                    return false;
                }
                let c2d2 = c2d2.unwrap();
                let mut pos2 = Curve2dEval::point_at(&c2d2, f2);
                if uclosed && ismodeu != 0 {
                    pos2.x += analysis::adjust_by_period(pos2.x, pos1.x, u_range);
                    if (pos2.x - pos1.x).abs() < PCONFUSION
                        && (found_u != 2 || pos1.x.abs() < uf.abs())
                    {
                        found_u = 2;
                        uf = pos1.x;
                    }
                }
                if vclosed && ismodeu == 0 {
                    pos2.y += analysis::adjust_by_period(pos2.y, pos1.y, v_range);
                    if (pos2.y - pos1.y).abs() < PCONFUSION
                        && (found_v != 2 || pos1.y.abs() < vf.abs())
                    {
                        found_v = 2;
                        vf = pos1.y;
                    }
                }
                i2 += 1;
            }
            i1 += 1;
        }

        // OCCT L2226-2234: pdn fixing RTS on offsets.
        if uf < su_f || uf > su_l {
            uf += analysis::adjust_to_period(uf, su_f, su_f + u_range);
        }
        if vf < sv_f || vf > sv_l {
            vf += analysis::adjust_to_period(vf, sv_f, sv_f + v_range);
        }

        // OCCT L2236-2242: create fictive grid and call ComposeShell to
        // insert a seam.
        let rts = Surface3::Trimmed(rcad_kernel::geom::TrimmedSurface::new(
            self.my_surf.as_ref().unwrap().surface().clone(),
            uf,
            uf + u_range,
            vf,
            vf + v_range,
        ));
        let grid: Vec<Vec<Surface3>> = vec![vec![rts.clone()]];
        let g = ShapeExtendCompositeSurface::new_with_grid_param(
            grid,
            crate::shhealing::shape_extend::status::ShapeExtendParametrisation::Natural,
        );
        let l_loc = 0u32;

        // OCCT L2245-2250: addition of non-manifold topology (the sequence is
        // walked without body — the aSeqNonManif shapes are re-added by the
        // ComposeShell result in OCCT).
        for _j in 1..=(a_seq_non_manif.len() as i32) {
            // (OCCT L2248-2249: the non-manifold shapes are handled by the
            // ComposeShell — kept for the call-shape parity.)
        }

        // OCCT L2252-2261.
        let mut comp_shell = ShapeFixComposeShellGap::new();
        comp_shell.init(g, l_loc, &tmp_f, CONFUSION); // myPrecision
        if self.base.my_context.is_none() {
            self.base.set_context(ShapeBuildReShape::new());
        }
        *comp_shell.closed_mode() = true;
        comp_shell.set_context(self.base.my_context.clone());
        comp_shell.set_max_tolerance(self.base.my_max_tol);
        comp_shell.perform();

        // OCCT L2263-2264: abv 03.07.00: CAX-IF TRJ4 — reset mySurf.
        self.my_surf = Some(ShapeAnalysisSurface::new(rts));

        // OCCT L2266.
        self.my_result = comp_shell.result();

        {
            let old = self.my_face.clone();
            let ctx = self.base.my_context.as_mut().unwrap();
            ctx.replace(brep, &old, &self.my_result);
        }

        // OCCT L2270-2300: remove small wires and/or faces generated by
        // ComposeShell (tests bugs step bug30052_4, de step_3 E6).
        let mut nb_faces = 0i32;
        for exp_f in topexp_explorer(brep, &self.my_result, ShapeType::Face) {
            let a_face = exp_f;
            let mut nb_wires = 0i32;
            for a_exp_w in topexp_explorer(brep, &a_face, ShapeType::Wire) {
                let mut a_sfw = ShapeFixWire::with_wire_face(brep, &a_exp_w, &a_face, self.base.my_precision);
                a_sfw.base.my_context = self.base.my_context.clone();
                if a_sfw.nb_edges() != 0 {
                    a_sfw.fix_small_all(brep, true, self.base.my_precision);
                }
                if a_sfw.nb_edges() == 0 {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.remove(brep, &a_exp_w);
                    }
                    continue;
                }
                nb_wires += 1;
            }
            if nb_wires == 0 {
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.remove(brep, &a_face);
                }
                continue;
            }
            nb_faces += 1;
        }

        {
            let ctx = self.base.my_context.as_mut().unwrap();
            self.my_result = ctx.apply(brep, &self.my_result, ShapeType::Shape);
        }
        for exp in topexp_explorer(brep, &self.my_result, ShapeType::Face) {
            if let Some(ctx) = self.base.my_context.as_mut() {
                self.my_face = ctx.apply(brep, &exp, ShapeType::Shape);
            }
            if self.my_face.is_null() {
                continue;
            }
            if nb_faces > 1 {
                self.fix_small_area_wire(brep, true);
                if let Some(ctx) = self.base.my_context.as_mut() {
                    let a_shape = ctx.apply(brep, &self.my_face, ShapeType::Shape);
                    if a_shape.is_null() {
                        continue;
                    }
                    self.my_face = a_shape;
                }
            }
            let my_face = self.my_face.clone();
            super::split_tool::brep_tools_update(brep, &my_face); //: p4
        }
        {
            let ctx = self.base.my_context.as_mut().unwrap();
            self.my_result = ctx.apply(brep, &self.my_result, ShapeType::Shape);
        }

        // OCCT L2324: Missing seam-edge added
        self.base
            .send_warning_own(&MessageMsg::from_key("FixAdvFace.FixMissingSeam.MSG0"));
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L2331-2394 — FixSmallAreaWire.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixSmallAreaWire (cxx L2331-2394): %14 pdn
    /// 24.02.99 PRO10109, USA60293 — fixes the wire on a face with small
    /// area.
    pub fn fix_small_area_wire(&mut self, brep: &mut BRep, the_is_remove_small_face: bool) -> bool {
        // OCCT L2333-2337.
        if self.base.my_context.is_some() {
            let f = self.my_face.clone();
            let a_shape = self.context_apply(brep, &f);
            self.my_face = a_shape;
        }

        let mut nb_removed = 0i32;
        let mut nb_wires = 0i32;

        // OCCT L2339-2344.
        let an_empty_copy = brep.empty_copied(&self.my_face);
        let mut a_face = an_empty_copy;
        a_face.orientation = Orientation::Forward;

        let a_tolerance3d = self.base.my_precision;
        for value in iter_subshapes(brep, &self.my_face, false, true) {
            let a_shape = value;
            if a_shape.shape_type() != ShapeType::Wire
                && a_shape.orientation != Orientation::Forward
                && a_shape.orientation != Orientation::Reversed
            {
                continue;
            }

            let a_wire = a_shape;
            let mut an_analyzer =
                ShapeAnalysisWire::new_from_wire(brep, &a_wire, &self.my_face, a_tolerance3d);
            if an_analyzer.check_small_area(brep, &a_wire) {
                // Null area wire detected, wire skipped
                self.base.send_warning(
                    &a_wire,
                    &MessageMsg::from_key("FixAdvFace.FixSmallAreaWire.MSG0"),
                );
                nb_removed += 1;
            } else {
                builder_add(brep, &a_face, &a_wire);
                nb_wires += 1;
            }
        }

        // OCCT L2372-2375.
        if nb_removed <= 0 {
            return false;
        }

        // OCCT L2377-2385.
        if nb_wires <= 0 {
            if the_is_remove_small_face && self.base.my_context.is_some() {
                let old = self.my_face.clone();
                let ctx = self.base.my_context.as_mut().unwrap();
                ctx.remove(brep, &old);
            }

            return false;
        }
        // OCCT L2386-2390.
        a_face.orientation = self.my_face.orientation;
        if let Some(ctx) = self.base.my_context.as_mut() {
            let old = self.my_face.clone();
            ctx.replace(brep, &old, &a_face);
        }

        self.my_face = a_face;
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L2478-2642 — FixLoopWire.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixLoopWire (cxx L2478-2642): detects if the wire
    /// has a loop and fixes the situation by splitting into a few parts; the
    /// ShapeExtend_DONE6 status is set when the wire had loops and was split.
    pub fn fix_loop_wire(&mut self, brep: &mut BRep, a_res_wires: &mut Vec<Shape>) -> bool {
        // OCCT L2480-2491.
        let mut a_map_vertices: Vec<(u64, u32)> = Vec::new();
        let mut a_map_vertex_edges: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        let mut a_map_small_edges: Vec<(u64, u32)> = Vec::new();
        let mut a_map_seem_edges: Vec<(u64, u32)> = Vec::new();
        if !self
            .my_fix_wire
            .my_analyzer
            .check_loop(brep, &mut a_map_vertices, &mut a_map_vertex_edges, &mut a_map_small_edges, &mut a_map_seem_edges)
        {
            return false;
        }

        // OCCT L2493-2494.
        let mut a_map_edges: Vec<(u64, u32)> = Vec::new();
        let mut a_seq_wires: Vec<Shape> = Vec::new();

        // OCCT L2496-2548: collecting wires from the common vertex belonging
        // to more than 2 edges.
        for i in 1..=(a_map_vertices.len() as i32) {
            let a_vert = a_map_vertices[(i - 1) as usize];
            let aledges = a_map_vertex_edges.get(&a_vert).cloned().unwrap_or_default();
            for liter in aledges.iter() {
                let edge = liter.clone();
                if a_map_edges.contains(&(edge.ptr_id(), edge.location)) {
                    continue;
                }

                let mut a_wire_data = WireData::new();
                a_wire_data.add_edge(&edge, 0);
                if a_map_seem_edges.contains(&(edge.ptr_id(), edge.location)) {
                    a_wire_data.add_edge(&shape_oriented(&edge, Orientation::Reversed), 0);
                }
                a_map_edges.push((edge.ptr_id(), edge.location));
                find_next(
                    &a_vert,
                    &edge,
                    &a_map_vertices,
                    &a_map_vertex_edges,
                    &a_map_small_edges,
                    &a_map_seem_edges,
                    &mut a_map_edges,
                    &mut a_wire_data,
                );
                if a_wire_data.nb_edges() == 1
                    && a_map_small_edges.contains(&(
                        a_wire_data.edge(1).ptr_id(),
                        a_wire_data.edge(1).location,
                    ))
                {
                    continue;
                }
                let a_wd_wire = a_wire_data.wire(brep);
                let (a_v1, a_v2) = topexp_vertices(brep, &a_wd_wire);

                if !a_v1.is_null()
                    && (a_v1.ptr_id(), a_v1.location) == (a_v2.ptr_id(), a_v2.location)
                {
                    let asfw_wire = a_wire_data.wire(brep);
                    let mut asfw = ShapeFixWire::new();
                    let asewd = WireData::new_from_wire(brep, &asfw_wire, true, false);
                    asfw.load_wire_data(brep, asewd);
                    asfw.fix_reorder(brep, false);
                    let awire2 = asfw.wire(brep);
                    a_res_wires.push(awire2);
                } else {
                    a_seq_wires.push(a_wire_data.wire(brep));
                }
            }
        }

        // OCCT L2550-2623.
        if a_seq_wires.len() == 1 {
            a_res_wires.push(a_seq_wires[0].clone());
        } else {
            // collecting whole wire from two not closed wires having two
            // common vertices.
            let mut i = 1usize;
            while i <= a_seq_wires.len() {
                let (a_v1, a_v2) = topexp_vertices(brep, &a_seq_wires[i - 1]);
                let a_wire = a_seq_wires[i - 1].clone();
                let mut j = i + 1;
                while j <= a_seq_wires.len() {
                    let (a_v21, a_v22) = topexp_vertices(brep, &a_seq_wires[(j - 1) as usize]);
                    let a_wire2 = a_seq_wires[(j - 1) as usize].clone();
                    if (a_v1.is_same(&a_v21) || a_v1.is_same(&a_v22))
                        && (a_v2.is_same(&a_v21) || a_v2.is_same(&a_v22))
                    {
                        let mut asewd = WireData::new_from_wire(brep, &a_wire, true, false);
                        asewd.add_wire(brep, &a_wire2, 0);
                        let mut asfw = ShapeFixWire::new();
                        asfw.load_wire_data(brep, asewd);
                        asfw.fix_reorder(brep, false);
                        a_res_wires.push(asfw.wire(brep));
                        a_seq_wires.remove(j - 1);
                        self.my_status |= encode_status(ShapeExtendStatus::Done7);
                        break;
                    }
                    j += 1;
                }
                if j <= a_seq_wires.len() {
                    a_seq_wires.remove(i - 1);
                    // OCCT L2583: aSeqWires.Remove(i--) — the outer index
                    // stays (the for-loop ++ is mirrored by the absence of
                    // the manual increment here).
                    continue;
                }
                i += 1;
            }
            if a_seq_wires.len() < 3 {
                for i in 1..=(a_seq_wires.len() as i32) {
                    a_res_wires.push(a_seq_wires[(i - 1) as usize].clone());
                }
            } else {
                // collecting wires having one common vertex
                for i in 1..=(a_seq_wires.len() as i32) {
                    let (mut a_v1, mut a_v2) = topexp_vertices(brep, &a_seq_wires[(i - 1) as usize]);
                    let mut a_wire = a_seq_wires[(i - 1) as usize].clone();
                    let mut j = (i + 1) as usize;
                    while j <= a_seq_wires.len() {
                        let (a_v21, a_v22) = topexp_vertices(brep, &a_seq_wires[(j - 1) as usize]);
                        let a_wire2 = a_seq_wires[(j - 1) as usize].clone();
                        if (a_v1.is_same(&a_v21) || a_v1.is_same(&a_v22))
                            || (a_v2.is_same(&a_v21) || a_v2.is_same(&a_v22))
                        {
                            let mut asewd = WireData::new_from_wire(brep, &a_wire, true, false);
                            asewd.add_wire(brep, &a_wire2, 0);
                            let mut asfw = ShapeFixWire::new();
                            asfw.load_wire_data(brep, asewd);
                            asfw.fix_reorder(brep, false);
                            a_wire = asfw.wire(brep);
                            let (nv1, nv2) = topexp_vertices(brep, &a_wire);
                            a_v1 = nv1;
                            a_v2 = nv2;
                            a_seq_wires.remove(j - 1);
                            self.my_status |= encode_status(ShapeExtendStatus::Done7);
                        } else {
                            j += 1;
                        }
                    }
                    a_res_wires.push(a_wire);
                }
            }
        }
        // OCCT L2624.
        let mut is_closed = true;

        // OCCT L2626-2639: checking that obtained wires is closed in 2D space.
        let is_not_plane = !matches!(
            self.my_surf.as_ref().map(|s| s.surface()),
            Some(Surface3::Plane(_))
        );
        if self.my_surf.is_some() && is_not_plane {
            let empty_copied = brep.empty_copied(&self.my_face);
            let mut tmp_face = empty_copied;
            tmp_face.orientation = Orientation::Forward;

            for i in 1..=(a_res_wires.len() as i32) {
                if !is_closed {
                    break;
                }
                let awire = a_res_wires[(i - 1) as usize].clone();
                is_closed = is_closed2d(brep, &tmp_face, &awire);
            }
        }

        // OCCT L2641.
        !a_res_wires.is_empty() && is_closed
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L2646-2729 — SplitEdge(sewd, num, param, vert,
    // preci, boxes).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::SplitEdge (cxx L2646-2729): splits the num-th
    /// edge of the wire by the parameter using the vertex.
    pub(crate) fn split_edge(
        &mut self,
        brep: &mut BRep,
        sewd: &mut WireData,
        num: i32,
        param: f64,
        vert: &Shape,
        preci: f64,
        boxes: &mut crate::shhealing::shape_fix::intersection_tool::ShapeBoxes2d,
    ) -> bool {
        // OCCT L2654-2657.
        let edge = sewd.edge(num);
        let mut new_e1 = Shape::null();
        let mut new_e2 = Shape::null();
        let a_tool = ShapeFixSplitTool::new();
        if !a_tool.split_edge(
            brep,
            &edge,
            param,
            vert,
            &self.my_face,
            &mut new_e1,
            &mut new_e2,
            preci,
            0.01 * preci,
        ) {
            return false;
        }

        // OCCT L2659-2671: change context.
        let mut wd = WireData::new();
        wd.add_edge(&new_e1, 0);
        wd.add_edge(&new_e2, 0);
        if let Some(ctx) = self.base.my_context.as_mut() {
            let w = wd.wire(brep);
            ctx.replace(brep, &edge, &w);
        }
        let w = wd.wire(brep);
        for e in topexp_explorer(brep, &w, ShapeType::Edge) {
            super::split_tool::brep_tools_update_edge(brep, &e);
        }

        // OCCT L2673-2682: change sewd and boxes.
        sewd.set_edge(&new_e1, num);
        if num == sewd.nb_edges() {
            sewd.add_edge(&new_e2, 0);
        } else {
            sewd.add_edge(&new_e2, num + 1);
        }

        boxes.remove(&(edge.ptr_id(), edge.location));
        rebind_edge_boxes(brep, &new_e1, &self.my_face, boxes);
        rebind_edge_boxes(brep, &new_e2, &self.my_face, boxes);
        // OCCT L2726.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L2733-2817 — SplitEdge(sewd, num, param1,
    // param2, vert, preci, boxes).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::SplitEdge (cxx L2733-2817): splits the num-th
    /// edge of the wire removing the segment (param1, param2).
    pub(crate) fn split_edge_two_params(
        &mut self,
        brep: &mut BRep,
        sewd: &mut WireData,
        num: i32,
        param1: f64,
        param2: f64,
        vert: &Shape,
        preci: f64,
        boxes: &mut crate::shhealing::shape_fix::intersection_tool::ShapeBoxes2d,
    ) -> bool {
        // OCCT L2742-2745.
        let edge = sewd.edge(num);
        let mut new_e1 = Shape::null();
        let mut new_e2 = Shape::null();
        let a_tool = ShapeFixSplitTool::new();
        if !a_tool.split_edge_two_params(
            brep,
            &edge,
            param1,
            param2,
            vert,
            &self.my_face,
            &mut new_e1,
            &mut new_e2,
            preci,
            0.01 * preci,
        ) {
            return false;
        }

        // OCCT L2747-2759: change context.
        let mut wd = WireData::new();
        wd.add_edge(&new_e1, 0);
        wd.add_edge(&new_e2, 0);
        if let Some(ctx) = self.base.my_context.as_mut() {
            let w = wd.wire(brep);
            ctx.replace(brep, &edge, &w);
        }
        let w = wd.wire(brep);
        for e in topexp_explorer(brep, &w, ShapeType::Edge) {
            super::split_tool::brep_tools_update_edge(brep, &e);
        }

        // OCCT L2761-2770: change sewd and boxes.
        sewd.set_edge(&new_e1, num);
        if num == sewd.nb_edges() {
            sewd.add_edge(&new_e2, 0);
        } else {
            sewd.add_edge(&new_e2, num + 1);
        }

        boxes.remove(&(edge.ptr_id(), edge.location));
        rebind_edge_boxes(brep, &new_e1, &self.my_face, boxes);
        rebind_edge_boxes(brep, &new_e2, &self.my_face, boxes);
        // OCCT L2814.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L2821-2825 — FixIntersectingWires.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixIntersectingWires (cxx L2821-2825).
    pub(crate) fn fix_intersecting_wires(&mut self, brep: &mut BRep) -> bool {
        let mut i_tool = ShapeFixIntersectionTool::new(
            self.base.my_context.clone(),
            self.base.my_precision,
            self.base.my_max_tol,
        );
        let mut face = self.my_face.clone();
        let res = i_tool.fix_intersecting_wires(brep, &mut face);
        self.my_face = face;
        res
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L2829-2901 — FixWiresTwoCoincEdges.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixWiresTwoCoincEdges (cxx L2829-2901): if a wire
    /// contains two coincident edges it must be removed; queries the status
    /// after Perform.
    pub(crate) fn fix_wires_two_coinc_edges(&mut self, brep: &mut BRep) -> bool {
        // OCCT L2831-2835.
        if self.base.my_context.is_some() {
            let f = self.my_face.clone();
            let s = self.context_apply(brep, &f);
            self.my_face = s;
        }

        // OCCT L2837-2842.
        let ori = self.my_face.orientation;
        let empty_copied = brep.empty_copied(&self.my_face);
        let mut face = empty_copied;
        face.orientation = Orientation::Forward;
        let mut nb_wires = 0i32;

        for value in iter_subshapes(brep, &self.my_face, false, true) {
            if value.shape_type() != ShapeType::Wire
                || (value.orientation != Orientation::Forward
                    && value.orientation != Orientation::Reversed)
            {
                continue;
            }
            nb_wires += 1;
        }
        if nb_wires < 2 {
            return false;
        }
        // OCCT L2858-2889.
        let mut is_fixed = false;
        for value in iter_subshapes(brep, &self.my_face, false, true) {
            if value.shape_type() != ShapeType::Wire
                || (value.orientation != Orientation::Forward
                    && value.orientation != Orientation::Reversed)
            {
                builder_add(brep, &face, &value);
                continue;
            }
            let wire = value;
            let sewd = WireData::new_from_wire(brep, &wire, true, false);
            if sewd.nb_edges() == 2 {
                let mut e1 = sewd.edge(1);
                let mut e2 = sewd.edge(2);
                e1.orientation = Orientation::Forward;
                e2.orientation = Orientation::Forward;
                if !shape_is_equal(&e1, &e2) {
                    builder_add(brep, &face, &wire);
                } else {
                    is_fixed = true;
                }
            } else {
                builder_add(brep, &face, &wire);
            }
        }
        // OCCT L2890-2898.
        if is_fixed {
            face.orientation = ori;
            if let Some(ctx) = self.base.my_context.as_mut() {
                let old = self.my_face.clone();
                ctx.replace(brep, &old, &face);
            }
            self.my_face = face;
        }

        is_fixed
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L2905-3010 — FixSplitFace.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixSplitFace (cxx L2905-3010): splits the face if
    /// there are more than one out wire, using the information after
    /// FixOrientation.
    pub(crate) fn fix_split_face(&mut self, brep: &mut BRep, map_wires: &WireListMap) -> bool {
        // OCCT L2909-2915.
        let mut faces: Vec<Shape> = Vec::new();
        let mut s = self.my_face.clone();
        if self.base.my_context.is_some() {
            let f = self.my_face.clone();
            s = self.context_apply(brep, &f);
        }
        let mut nb_wires = 0i32;
        let mut nb_wires_new = 0i32;
        for value in iter_subshapes(brep, &s, false, true) {
            let a_shape = value;
            if a_shape.shape_type() != ShapeType::Wire
                || (a_shape.orientation != Orientation::Forward
                    && a_shape.orientation != Orientation::Reversed)
            {
                continue;
            }
            let wire = a_shape;
            nb_wires += 1;
            if !map_wires.contains_key(&(wire.ptr_id(), wire.location)) {
                continue;
            }
            // if wire not closed --> stop split and return false
            // OCCT L2930-2935.
            let sewd = WireData::new_from_wire(brep, &wire, true, false);
            let nb_edges = sewd.nb_edges();
            if nb_edges == 0 {
                continue;
            }
            //
            // OCCT L2937-2946.
            let e1 = sewd.edge(1);
            let e2 = sewd.edge(nb_edges);
            let sae = ShapeAnalysisEdge::new();
            let v1 = sae.first_vertex(brep, &e1);
            let v2 = sae.last_vertex(brep, &e2);
            if !v1.is_same(&v2) {
                return false;
            }
            // OCCT L2947-2952: create face.
            let empty_copied = brep.empty_copied(&s);
            let mut tmp_face = empty_copied;
            tmp_face.orientation = Orientation::Forward;
            builder_add(brep, &tmp_face, &wire);
            nb_wires_new += 1;
            let int_wires = map_wires.get(&(wire.ptr_id(), wire.location)).cloned().unwrap_or_default();
            for liter in int_wires.iter() {
                // OCCT L2957-2961.
                let a_shape_empty_copied = brep.empty_copied(&tmp_face);
                let mut a_face = a_shape_empty_copied;
                a_face.orientation = Orientation::Forward;
                builder_add(brep, &a_face, liter);
                let clas = FClass2dTopol::new(Arc::new(brep.clone()), &a_face, PCONFUSION);
                let staout = clas.perform_infinite_point();
                if staout == rcad_kernel::topods::State::In {
                    builder_add(brep, &tmp_face, liter);
                } else {
                    builder_add(brep, &tmp_face, &shape_oriented(liter, Orientation::Reversed));
                }
                nb_wires_new += 1;
            }
            if !self.my_fwd {
                tmp_face.orientation = Orientation::Reversed;
            }
            faces.push(tmp_face);
        }

        // OCCT L2981-2984.
        if nb_wires != nb_wires_new {
            return false;
        }

        // OCCT L2986-3007.
        if faces.len() > 1 {
            let comp = brep.add_tcompound(Vec::new());
            for i in 1..=(faces.len() as i32) {
                builder_add(brep, &comp, &faces[(i - 1) as usize]);
            }
            self.my_result = comp;

            if let Some(ctx) = self.base.my_context.as_mut() {
                let old = self.my_face.clone();
                ctx.replace(brep, &old, &self.my_result);
            }

            for exp in topexp_explorer(brep, &self.my_result, ShapeType::Face) {
                self.my_face = exp;
                let my_face = self.my_face.clone();
                super::split_tool::brep_tools_update(brep, &my_face);
            }
            return true;
        }

        // OCCT L3009.
        false
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L3101-3259 — FixPeriodicDegenerated.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixPeriodicDegenerated (cxx L3101-3259): fixes
    /// the topology for the case when the face is composed by a single wire
    /// belting a periodic surface — a degenerated edge is reconstructed in
    /// the degenerated pole of the surface.
    pub(crate) fn fix_periodic_degenerated(&mut self, brep: &mut BRep) -> bool {
        // OCCT L3103-3111: prepare fix routine.
        if self.base.my_context.is_some() {
            let f = self.my_face.clone();
            let a_sh = self.context_apply(brep, &f);
            self.my_face = a_sh;
        }

        // OCCT L3113-3129: check if fix can be applied on the passed face.
        let mut a_wire_seq: Vec<Shape> = Vec::new();
        for value in iter_subshapes(brep, &self.my_face, false, true) {
            let a_sub_sh = value;
            if a_sub_sh.shape_type() != ShapeType::Wire
                || (a_sub_sh.orientation != Orientation::Forward
                    && a_sub_sh.orientation != Orientation::Reversed)
            {
                continue;
            }
            a_wire_seq.push(a_sub_sh);
        }

        // OCCT L3131-3140.
        let a_nb_wires = a_wire_seq.len() as i32;
        let a_surface = super::split_tool::brep_tool_surface(brep, &self.my_face);

        // Only single wires on conical surfaces are checked.
        let is_cone = matches!(a_surface, Some(Surface3::Cone(_)));
        if a_nb_wires != 1 || a_surface.is_none() || !is_cone {
            return false;
        }

        // OCCT L3142-3162.
        let mut a_sole_wire = a_wire_seq[0].clone();
        let (mut a_min_loop_u, mut a_max_loop_u, mut a_min_loop_v, mut a_max_loop_v) =
            (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut is_u_decrease = false;

        let a_cone_surf = match &a_surface {
            Some(Surface3::Cone(c)) => c.clone(),
            _ => return false,
        };
        let is_conic_loop = is_periodic_conical_loop(
            &a_cone_surf,
            brep,
            &a_sole_wire,
            self.base.my_precision,
            &mut a_min_loop_u,
            &mut a_max_loop_u,
            &mut a_min_loop_v,
            &mut a_max_loop_v,
            &mut is_u_decrease,
        );

        if !is_conic_loop {
            return false;
        }

        // OCCT L3164-3183: retrieve apex.
        // Get base circle of the conical surface (the circle it was built
        // from): aConeBaseCrv = aConeSurf->VIso(0.0).
        let mut sas = ShapeAnalysisSurface::new(Surface3::Cone(a_cone_surf.clone()));
        let a_cone_base_crv = sas.viso(0.0);
        let a_cone_base_r = match &a_cone_base_crv {
            Some(Curve3::Circle(c)) => c.radius,
            // OCCT down_cast null — the arm cannot trigger for a cone.
            _ => a_cone_surf.radius,
        };

        // Retrieve conical props.
        let a_semi_angle = a_cone_surf.half_angle_rad;

        if a_semi_angle.abs() <= CONFUSION {
            return false; // Bad surface
        }

        // OCCT L3182-3183: find the V parameter of the apex.
        let a_cone_base_h = a_cone_base_r / a_semi_angle.sin();
        let an_apex_v = -a_cone_base_h;

        // OCCT L3185-3186: get apex vertex (BRepBuilderAPI_MakeVertex).
        let apex_pnt = a_cone_surf.apex - a_cone_surf.axis * (a_cone_surf.radius / a_semi_angle_rad_tan(a_semi_angle));
        let mut b = BRepBuilder::new();
        let an_apex = b.add_vertex(brep, apex_pnt, CONFUSION);

        // OCCT L3192-3202.
        // Check the positional relationship between the initial wire and the
        // apex line in 2D.
        if (an_apex_v - a_min_loop_v).abs() <= self.base.my_precision
            || (an_apex_v - a_max_loop_v).abs() <= self.base.my_precision
            || (an_apex_v < a_max_loop_v && an_apex_v > a_min_loop_v)
        {
            return false;
        }

        // OCCT L3204-3224.
        let an_apex_curve2d;
        // Apex curve below the wire
        if an_apex_v < a_min_loop_v {
            an_apex_curve2d = Curve2d::Line(Line2d {
                origin: DVec2::new(a_min_loop_u, an_apex_v),
                direction: DVec2::X,
            });
            if !is_u_decrease {
                a_sole_wire.orientation = match a_sole_wire.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };
            }
        // Apex curve above the wire
        } else if an_apex_v > a_max_loop_v {
            an_apex_curve2d = Curve2d::Line(Line2d {
                origin: DVec2::new(a_max_loop_u, an_apex_v),
                direction: -DVec2::X,
            });
            if is_u_decrease {
                a_sole_wire.orientation = match a_sole_wire.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };
            }
        } else {
            return false;
        }

        // OCCT L3226-3232: create degenerated edge & wire for apex.
        let an_apex_fwd = shape_oriented(&an_apex, Orientation::Forward);
        let an_apex_rev = shape_oriented(&an_apex, Orientation::Reversed);
        let mut an_apex_edge = b.add_edge(brep, None, an_apex_fwd, an_apex_rev, [0.0, 0.0]);
        b.update_edge_pcurve(brep, an_apex_edge.clone(), an_apex_curve2d, self.my_face.clone(), self.base.my_precision);
        set_edge_degenerated(brep, &an_apex_edge, true);
        b.set_edge_range(brep, an_apex_edge.clone(), 0.0, (a_max_loop_u - a_min_loop_u).abs());
        let an_apex_wire = brep.add_twire(vec![an_apex_edge]);

        // OCCT L3234-3256.
        let mut a_new_wire_seq: Vec<Shape> = Vec::new();
        a_new_wire_seq.push(a_sole_wire.clone());
        a_new_wire_seq.push(an_apex_wire);

        // Assemble new face.
        let a_new_face_empty = brep.empty_copied(&self.my_face);
        let mut a_new_face = a_new_face_empty;
        a_new_face.orientation = Orientation::Forward;
        let a_face_builder = BRepBuilder::new();
        let _ = a_face_builder;
        for i in 1..=(a_new_wire_seq.len() as i32) {
            let a_new_wire = a_new_wire_seq[(i - 1) as usize].clone();
            builder_add(brep, &a_new_face, &a_new_wire);
        }
        a_new_face.orientation = self.my_face.orientation;

        // Adjust the resulting state of the healing tool.
        self.my_result = a_new_face;
        if let Some(ctx) = self.base.my_context.as_mut() {
            let old = self.my_face.clone();
            let res = self.my_result.clone();
            ctx.replace(brep, &old, &res);
        }

        // OCCT L3258.
        true
    }
}

fn a_semi_angle_rad_tan(a_semi_angle: f64) -> f64 {
    a_semi_angle.tan()
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L1652-1718 — CheckWire (static).
// ---------------------------------------------------------------------------

/// OCCT static CheckWire (cxx L1652-1718): tests whether the wire is opened
/// on the period of the periodical surface (:i7 abv 18 Sep 98).
#[allow(clippy::too_many_arguments)]
pub(crate) fn check_wire(
    brep: &mut BRep,
    wire: &Shape,
    face: &Shape,
    d_u: f64,
    d_v: f64,
    isuopen: &mut i32,
    isvopen: &mut i32,
    isdeg: &mut bool,
) -> bool {
    let mut vec = DVec2::ZERO;
    let sae = ShapeAnalysisEdge::new();

    *isuopen = 0;
    *isvopen = 0;
    *isdeg = true;
    for value in iter_subshapes(brep, wire, true, true) {
        let edge = value;
        if !brep_tool_degenerated(&edge) {
            *isdeg = false;
        }
        let mut c2d: Option<Curve2d> = None;
        let mut f = 0.0f64;
        let mut l = 0.0f64;
        if !sae.pcurve_face(brep, &edge, face, &mut c2d, &mut f, &mut l, true) {
            return false;
        }
        let c2d = c2d.unwrap();
        let pl = Curve2dEval::point_at(&c2d, l);
        let pf = Curve2dEval::point_at(&c2d, f);
        vec += pl - pf;
    }

    // OCCT L1683-1698.
    let a_delta = vec.x.abs() - d_u;
    if a_delta.abs() < 0.1 * d_u {
        if vec.x > 0.0 {
            *isuopen = 1;
        } else {
            *isuopen = -1;
        }
    } else {
        *isuopen = 0;
    }

    // OCCT L1700-1715.
    let a_delta = vec.y.abs() - d_v;
    if a_delta.abs() < 0.1 * d_v {
        if vec.y > 0.0 {
            *isvopen = 1;
        } else {
            *isvopen = -1;
        }
    } else {
        *isvopen = 0;
    }

    // OCCT L1717.
    *isuopen != 0 || *isvopen != 0
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L2398-2456 — FindNext (static).
// ---------------------------------------------------------------------------

/// OCCT static FindNext (cxx L2398-2456): walks the vertex-edge adjacency to
/// collect the wire pieces through the loop vertices.
#[allow(clippy::too_many_arguments)]
pub(crate) fn find_next(
    a_vert: &(u64, u32),
    ainit_edge: &Shape,
    a_map_vertices: &[(u64, u32)],
    a_map_vertex_edges: &HashMap<(u64, u32), Vec<Shape>>,
    a_map_small_edges: &[(u64, u32)],
    a_map_seem_edges: &[(u64, u32)],
    a_map_edges: &mut Vec<(u64, u32)>,
    a_wire_data: &mut WireData,
) {
    // OCCT L2409-2419: the other vertex of the edge.
    let mut anext_vert = *a_vert;
    let mut is_find = false;
    if let TShape::Edge(ed) = ainit_edge.data.as_ref() {
        for v in [&ed.first, &ed.last] {
            let mut vv = v.clone();
            if ainit_edge.orientation == Orientation::Reversed {
                vv.orientation = match vv.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };
            }
            let key = (vv.ptr_id(), vv.location);
            if key != *a_vert {
                is_find = true;
                anext_vert = key;
            }
        }
    }

    // OCCT L2421-2428.
    if !is_find && !a_map_small_edges.contains(&(ainit_edge.ptr_id(), ainit_edge.location)) {
        return;
    }
    if is_find && a_map_vertices.contains(&anext_vert) {
        return;
    }

    // OCCT L2430-2455.
    let aledges = a_map_vertex_edges.get(&anext_vert).cloned().unwrap_or_default();
    let mut is_find2 = false;
    for liter in aledges.iter() {
        if is_find2 {
            break;
        }
        let lkey = (liter.ptr_id(), liter.location);
        if !a_map_edges.contains(&lkey) && !occt_is_same(liter, ainit_edge) {
            let anext_edge = liter.clone();
            a_wire_data.add_edge(&anext_edge, 0);
            if a_map_seem_edges.contains(&lkey) {
                a_wire_data.add_edge(&shape_oriented(&anext_edge, Orientation::Reversed), 0);
            }
            is_find2 = true;
            a_map_edges.push(lkey);
            find_next(
                &anext_vert,
                &anext_edge,
                a_map_vertices,
                a_map_vertex_edges,
                a_map_small_edges,
                a_map_seem_edges,
                a_map_edges,
                a_wire_data,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L2458-2474 — isClosed2D (static).
// ---------------------------------------------------------------------------

/// OCCT static isClosed2D (cxx L2458-2474).
pub(crate) fn is_closed2d(brep: &mut BRep, a_face: &Shape, a_wire: &Shape) -> bool {
    let mut is_closed = true;
    let mut asaw = ShapeAnalysisWire::new_from_wire(brep, a_wire, a_face, CONFUSION);
    let mut i = 1i32;
    while i <= asaw.nb_edges() && is_closed {
        let edge1 = asaw.wire_data().unwrap().edge(i);
        // checking that wire is closed in 2D space with tolerance of vertex.
        let sae = ShapeAnalysisEdge::new();
        let v1 = sae.first_vertex(brep, &edge1);
        asaw.set_precision(brep_tool_tolerance(&v1));
        asaw.check_gap2d(brep, i);
        is_closed = asaw.last_check_status(ShapeExtendStatus::Ok);
        i += 1;
    }
    is_closed
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L3018-3097 — IsPeriodicConicalLoop (static).
// ---------------------------------------------------------------------------

/// OCCT static IsPeriodicConicalLoop (cxx L3018-3097): checks whether the
/// passed wire makes up a periodic loop on the passed conical surface.
#[allow(clippy::too_many_arguments)]
pub(crate) fn is_periodic_conical_loop(
    the_surf: &rcad_kernel::geom::ConicalSurface,
    brep: &mut BRep,
    the_wire: &Shape,
    the_tolerance: f64,
    the_min_u: &mut f64,
    the_max_u: &mut f64,
    the_min_v: &mut f64,
    the_max_v: &mut f64,
    is_u_decrease: &mut bool,
) -> bool {
    let a_sae = ShapeAnalysisEdge::new();

    let mut a_cumul_delta_u = 0.0f64;
    let mut a_cumul_delta_u_abs = 0.0f64;
    let mut a_min_u = f64::MAX;
    let mut a_min_v = a_min_u;
    let mut a_max_u = -a_min_u;
    let mut a_max_v = a_max_u;

    // OCCT L3042-3085: iterate over the edges.
    let _ = brep;
    let _ = the_wire;
    for value in iter_subshapes(brep, the_wire, false, true) {
        let a_current_edge = value;
        let mut a_c2d: Option<Curve2d> = None;
        let mut a_pfirst = 0.0f64;
        let mut a_plast = 0.0f64;

        // OCCT L3049: aSAE.PCurve(aCurrentEdge, theSurf, aLoc, aC2d,
        // aPFirst, aPLast, true) — the (surface, location) overload; the
        // location is the default (identity) TopLoc_Location.
        a_sae.pcurve_surface(
            brep,
            &a_current_edge,
            &Surface3::Cone(the_surf.clone()),
            0,
            &mut a_c2d,
            &mut a_pfirst,
            &mut a_plast,
            true,
        );
        let a_c2d = match a_c2d {
            Some(c) => c,
            // OCCT L3051-3054.
            None => return false,
        };

        let a_uvfirst = Curve2dEval::point_at(&a_c2d, a_pfirst);
        let a_uvlast = Curve2dEval::point_at(&a_c2d, a_plast);

        let (a_ufirst, a_ulast) = (a_uvfirst.x, a_uvlast.x);
        let (a_vfirst, a_vlast) = (a_uvfirst.y, a_uvlast.y);

        let a_cur_max_u = a_ufirst.max(a_ulast);
        let a_cur_min_u = a_ufirst.min(a_ulast);
        let a_cur_max_v = a_vfirst.max(a_vlast);
        let a_cur_min_v = a_vfirst.min(a_vlast);

        if a_cur_min_u < a_min_u {
            a_min_u = a_cur_min_u;
        }
        if a_cur_max_u > a_max_u {
            a_max_u = a_cur_max_u;
        }
        if a_cur_min_v < a_min_v {
            a_min_v = a_cur_min_v;
        }
        if a_cur_max_v > a_max_v {
            a_max_v = a_cur_max_v;
        }

        let a_delta_u = a_ulast - a_ufirst;

        a_cumul_delta_u += a_delta_u;
        a_cumul_delta_u_abs += a_delta_u.abs();
    }

    *the_min_u = a_min_u;
    *the_max_u = a_max_u;
    *the_min_v = a_min_v;
    *the_max_v = a_max_v;
    *is_u_decrease = a_cumul_delta_u < 0.0;

    // OCCT L3093-3096.
    let is_2pi_delta = (a_cumul_delta_u_abs - 2.0 * std::f64::consts::PI).abs() <= the_tolerance;
    let is_around_apex = (*the_max_u - *the_min_u).abs() > 2.0 * std::f64::consts::PI - the_tolerance;

    is_2pi_delta && is_around_apex
}


// ---------------------------------------------------------------------------
// The box-rebind block shared by SplitEdge x2 (cxx L2684-2725 /
// L2772-2813).
// ---------------------------------------------------------------------------

/// OCCT cxx L2685-2706 (repeated for newE2): the pcurve box of the new edge
/// bound in the boxes map.
pub(crate) fn rebind_edge_boxes(
    brep: &mut BRep,
    new_e: &Shape,
    face: &Shape,
    boxes: &mut crate::shhealing::shape_fix::intersection_tool::ShapeBoxes2d,
) {
    let (s, l) = brep_tool_surface_loc(brep, face);
    let s = match s {
        Some(s) => s,
        None => return,
    };
    let sae = ShapeAnalysisEdge::new();
    let mut c2d: Option<Curve2d> = None;
    let mut cf = 0.0f64;
    let mut cl = 0.0f64;
    if sae.pcurve_surface(brep, new_e, &s, l, &mut c2d, &mut cf, &mut cl, false) {
        let mut box2d = rcad_kernel::math::bnd::BndBox2d::new();
        let c = c2d.as_ref().unwrap();
        let domain = c.default_domain();
        let (a_first, a_last) = (domain[0], domain[1]);
        if matches!(c, Curve2d::BSpline(_)) && (cf < a_first || cl > a_last) {
            // pdn avoiding problems with segment in Bnd_Box
            bnd_lib_add2d_curve(brep, c, a_first, a_last, &mut box2d);
        } else {
            bnd_lib_add2d_curve(brep, c, cf, cl, &mut box2d);
        }
        boxes.insert((new_e.ptr_id(), new_e.location), box2d);
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_ComposeShell — W3 docket row, C-on-demand (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `ShapeFix_ComposeShell` (ShapeFix_ComposeShell.cxx, 3,606 LOC; the
/// docket C-on-demand row — zero oracle in the 36-case heal face) — GAP
/// carrier reduced to the call shape of `FixMissingSeam` (cxx L2252-2266):
/// Init / ClosedMode / SetContext / SetMaxTolerance / Perform / Result.
/// `Perform` keeps OCCT's "nothing composed" path: the result is the input
/// face unchanged.  Replaced wholesale when the ComposeShell tranche lands.
pub struct ShapeFixComposeShellGap {
    my_face: Shape,
    my_closed_mode: bool,
}

impl ShapeFixComposeShellGap {
    /// OCCT ShapeFix_ComposeShell() (the default constructor).
    pub fn new() -> Self {
        ShapeFixComposeShellGap {
            my_face: Shape::null(),
            my_closed_mode: false,
        }
    }

    /// OCCT ShapeFix_ComposeShell::Init(Grid, L, T, prec) (cxx L214-246) —
    /// GAP: the arguments keep the call shape; only the shape is stored.
    pub fn init(
        &mut self,
        _grid: ShapeExtendCompositeSurface,
        _loc: u32,
        t: &Shape,
        _preci: f64,
    ) {
        self.my_face = t.clone();
    }

    /// OCCT ShapeFix_ComposeShell::ClosedMode() (lxx).
    pub fn closed_mode(&mut self) -> &mut bool {
        &mut self.my_closed_mode
    }

    /// OCCT ShapeFix_ComposeShell::SetContext (cxx L97-100).
    pub fn set_context(&mut self, _context: Option<ShapeBuildReShape>) {}

    /// OCCT ShapeFix_Root::SetMaxTolerance.
    pub fn set_max_tolerance(&mut self, _maxtol: f64) {}

    /// OCCT ShapeFix_ComposeShell::Perform (cxx L248-...) — GAP: the
    /// "nothing composed" path.
    pub fn perform(&mut self) -> bool {
        false
    }

    /// OCCT ShapeFix_ComposeShell::Result (lxx) — the no-result path of the
    /// GAP: the input face passes through unchanged.
    pub fn result(&self) -> Shape {
        self.my_face.clone()
    }
}

impl Default for ShapeFixComposeShellGap {
    fn default() -> Self {
        Self::new()
    }
}

