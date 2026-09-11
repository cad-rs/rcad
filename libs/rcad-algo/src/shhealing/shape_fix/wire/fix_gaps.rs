//! OCCT `ShapeFix_Wire` — the gap-fixing methods of `ShapeFix_Wire_1.cxx`
//! (L1-2006): `FixGaps3d` / `FixGaps2d` / `FixGap3d` / `FixGap2d` (the szv
//! 19.08.99 methods for fixing gaps between edges — 3d curves and pcurves).
//! The impl submodule of [`super::ShapeFixWire`] (the OCCT continued-file
//! convention); the `AdjustOnPeriodic3d`/`AdjustOnPeriodic2d` statics live
//! in `wire_statics.rs`.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, OffsetCurve2d, OffsetCurve3, TrimmedCurve2, TrimmedCurve3};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType, TShape};

use crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::transfer_parameters_proj::ShapeAnalysisTransferParametersProj;
use crate::shhealing::shape_build::edge::{builder_range_on_face, ShapeBuildEdge};
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::shape_tolerance::ShapeFixShapeTolerance;
use crate::shhealing::shape_fix::wire::wire_statics::{
    adjust_on_periodic2d, adjust_on_periodic3d, builder_update_edge_curve3d_value,
    pcurve_range_on_face,
};
use crate::shhealing::shape_fix::wire::ShapeFixWire;


impl ShapeFixWire {
    // OCCT ShapeFix_Wire_1.cxx L76-97 — FixGaps3d.
    /// OCCT ShapeFix_Wire::FixGaps3d (Wire_1.cxx L76-97): fixes gaps between
    /// the ends of the 3d curves on adjacent edges; myPrecision is used to
    /// detect the gaps.
    pub fn fix_gaps3d(&mut self, brep: &mut BRep) -> bool {
        self.my_status_gaps3d = encode_status(ShapeExtendStatus::Ok);
        // if ( !IsReady() ) return false;

        let start = if self.my_closed_mode { 1 } else { 2 };
        if self.my_fix_gaps_by_ranges {
            for i in start..=self.nb_edges() {
                self.fix_gap3d(brep, i, false);
                self.my_status_gaps3d |= self.my_last_fix_status;
            }
        }
        for i in start..=self.nb_edges() {
            self.fix_gap3d(brep, i, true);
            self.my_status_gaps3d |= self.my_last_fix_status;
        }

        self.status_gaps3d(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire_1.cxx L101-122 — FixGaps2d.
    /// OCCT ShapeFix_Wire::FixGaps2d (Wire_1.cxx L101-122): fixes gaps
    /// between the ends of the pcurves on adjacent edges.
    pub fn fix_gaps2d(&mut self, brep: &mut BRep) -> bool {
        self.my_status_gaps2d = encode_status(ShapeExtendStatus::Ok);
        //  if ( !IsReady() ) return false;

        let start = if self.my_closed_mode { 1 } else { 2 };
        if self.my_fix_gaps_by_ranges {
            for i in start..=self.nb_edges() {
                self.fix_gap2d(brep, i, false);
                self.my_status_gaps2d |= self.my_last_fix_status;
            }
        }
        for i in start..=self.nb_edges() {
            self.fix_gap2d(brep, i, true);
            self.my_status_gaps2d |= self.my_last_fix_status;
        }

        self.status_gaps2d(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire_1.cxx L154-861 — FixGap3d.
    /// OCCT ShapeFix_Wire::FixGap3d(num, convert) (Wire_1.cxx L154-861):
    /// fixes the gap between the ends of the 3d curves on the (num-1)-th and
    /// num-th edges.
    pub fn fix_gap3d(&mut self, brep: &mut BRep, num: i32, convert: bool) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        //  if ( !IsReady() ) return false;

        //=============
        // First phase: analysis whether the problem (gap) exists
        //=============

        if self.base.my_context.is_none() {
            self.base.set_context(crate::shhealing::shape_build::reshape::ShapeBuildReShape::new());
        }

        let preci = self.base.my_precision;

        let sbwd_nb = self.my_analyzer.wire_data().unwrap().nb_edges();
        let n2 = if num > 0 { num } else { sbwd_nb };
        let n1 = if n2 > 1 { n2 - 1 } else { sbwd_nb };
        // smh#8
        let raw1 = self.my_analyzer.wire_data().unwrap().edge(n1);
        let raw2 = self.my_analyzer.wire_data().unwrap().edge(n2);
        let e1 = self.context_apply(brep, &raw1);
        let e2 = self.context_apply(brep, &raw2);

        // Retrieve curves on edges
        let sae = ShapeAnalysisEdge::new();
        let mut c1: Option<Curve3> = None;
        let mut c2: Option<Curve3> = None;
        let mut cfirst1 = 0.0f64;
        let mut clast1 = 0.0f64;
        let mut cfirst2 = 0.0f64;
        let mut clast2 = 0.0f64;
        if !sae.curve3d(brep, &e1, &mut c1, &mut cfirst1, &mut clast1, false)
            || !sae.curve3d(brep, &e2, &mut c2, &mut cfirst2, &mut clast2, false)
        {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let c1 = c1.unwrap();
        let c2 = c2.unwrap();

        // Check gap in 3d space
        let cpnt1 = CurveEval::point_at(&c1, clast1);
        let cpnt2 = CurveEval::point_at(&c2, cfirst2);
        let gap = cpnt1.distance(cpnt2);
        if !convert && gap <= preci {
            return false;
        }

        //=============
        // Second phase: collecting data necessary for further analysis
        //=============

        let reversed1 = e1.orientation == Orientation::Reversed;
        let reversed2 = e2.orientation == Orientation::Reversed;

        let v1 = sae.last_vertex(brep, &e1);
        let v2 = sae.first_vertex(brep, &e2);
        let mut vpnt = if v1.is_same(&v2) {
            brep_tool_pnt_of(&v1)
        } else {
            (brep_tool_pnt_of(&v1) + brep_tool_pnt_of(&v2)) * 0.5
        };

        let mut first1;
        let mut last1;
        let mut first2;
        let mut last2;
        if reversed1 {
            first1 = clast1;
            last1 = cfirst1;
        } else {
            first1 = cfirst1;
            last1 = clast1;
        }
        if reversed2 {
            first2 = clast2;
            last2 = cfirst2;
        } else {
            first2 = cfirst2;
            last2 = clast2;
        }

        let mut c1v = c1.clone();
        let mut c2v = c2.clone();

        // Extract basic curves from trimmed and offset
        let mut basic = false;
        let mut trimmed1 = false;
        let mut offset1 = false;
        let mut offval1 = DVec3::ZERO;
        while !basic {
            let next = c1v.clone();
            match next {
                Curve3::Trimmed(trm) => {
                    c1v = (*trm.curve).clone();
                    trimmed1 = true;
                }
                Curve3::Offset(oc) => {
                    c1v = (*oc.basis).clone();
                    offval1 += oc.offset_dir * oc.offset_distance;
                    offset1 = true;
                }
                _ => basic = true,
            }
        }
        basic = false;
        let mut trimmed2 = false;
        let mut offset2 = false;
        let mut offval2 = DVec3::ZERO;
        while !basic {
            let next = c2v.clone();
            match next {
                Curve3::Trimmed(trm) => {
                    c2v = (*trm.curve).clone();
                    trimmed2 = true;
                }
                Curve3::Offset(oc) => {
                    c2v = (*oc.basis).clone();
                    offval2 += oc.offset_dir * oc.offset_distance;
                    offset2 = true;
                }
                _ => basic = true,
            }
        }
        // Restore offset curves
        if offset1 {
            c1v = Curve3::Offset(OffsetCurve3 {
                basis: Box::new(c1v),
                offset_distance: offval1.length(),
                offset_dir: offval1.normalize_or_zero(),
            });
        }
        if offset2 {
            c2v = Curve3::Offset(OffsetCurve3 {
                basis: Box::new(c2v),
                offset_distance: offval2.length(),
                offset_dir: offval2.normalize_or_zero(),
            });
        }

        let mut done1 = false;
        let mut done2 = false;

        if convert {
            // Check that the gap satisfies the precision — in this case no
            // conversion is produced
            if cpnt1.distance(vpnt) < preci && cpnt2.distance(vpnt) < preci {
                return false;
            }

            let mut bsp1: Option<rcad_kernel::geom::BSplineCurve3> = None;
            let mut bsp2: Option<rcad_kernel::geom::BSplineCurve3> = None;

            // iterate on curves
            let nbcurv = if n1 == n2 { 1 } else { 2 };
            for j in 1..=nbcurv {
                // bool trim = false;  // skl
                let c: &Curve3;
                let first: f64;
                let last: f64;
                if j == 1 {
                    if cpnt1.distance(vpnt) < preci {
                        if n1 == n2 && cpnt2.distance(vpnt) < preci {
                            continue;
                        }
                        if n1 != n2 {
                            continue;
                        }
                    }
                    c = &c1v;
                    first = first1;
                    last = last1; /*trim = trimmed1;*/ // skl
                } else {
                    if cpnt2.distance(vpnt) < preci {
                        continue;
                    }
                    c = &c2v;
                    first = first2;
                    last = last2; /*trim = trimmed2;*/ // skl
                }

                let mut bsp: Option<rcad_kernel::geom::BSplineCurve3> = None;

                // Convert the curve to a bspline
                match c {
                    Curve3::BSpline(bs0) => {
                        let mut bs = bs0.clone();
                        // take segment if trim and range differ
                        let mut fbsp = CurveEval::default_domain(&Curve3::BSpline(bs.clone()))[0];
                        let mut lbsp = CurveEval::default_domain(&Curve3::BSpline(bs.clone()))[1];
                        let mut segment = false;
                        if first > fbsp {
                            fbsp = first;
                            segment = true;
                        }
                        if last < lbsp {
                            lbsp = last;
                            segment = true;
                        }
                        if segment {
                            // OCCT: GeomConvert::SplitBSplineCurve(bsp, fbsp,
                            // lbsp, Confusion) — the TKGeomBase split GAP
                            // keeps the no-result path.
                            bsp = geom_convert_split_bspline(&bs, fbsp, lbsp, CONFUSION);
                        } else {
                            bsp = Some(bs);
                        }
                    }
                    Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::Hyperbola(_) | Curve3::Parabola(_) => {
                        // OCCT L361-372: Approx_Curve3d(Conv, C1, 9, 1000) —
                        // the TKGeomBase approximation GAP keeps the
                        // `!IsDone()` path.
                        bsp = None;
                    }
                    _ => {
                        // Restore trim for pcurve
                        // OCCT L376-402: the try/catch arm (bridge #5).
                        let sac = ShapeAnalysisCurve;
                        let (tf, tl) = if !sac.is_periodic(c) {
                            (
                                first.max(CurveEval::default_domain(c)[0]),
                                last.min(CurveEval::default_domain(c)[1]),
                            )
                        } else {
                            (first, last)
                        };
                        // OCCT: GeomConvert::CurveToBSplineCurve(tc) — the
                        // TKGeomBase conversion GAP keeps the catch path.
                        bsp = geom_convert_curve_to_bspline(
                            &Curve3::Trimmed(TrimmedCurve3::new(c.clone(), tf, tl)),
                        );
                    }
                }

                if j == 1 {
                    bsp1 = bsp;
                } else {
                    bsp2 = bsp;
                }
            }

            // Take the curve ends if the conversion failed
            if bsp1.is_none() {
                vpnt = cpnt1;
            } else if bsp2.is_none() {
                vpnt = cpnt2;
            }

            if let Some(bs) = bsp1.as_mut() {
                if bs.degree == 1 {
                    bs.increase_degree(2); // gka
                }
                if n1 == n2 {
                    set_pole3(bs, 1, vpnt);
                    let np = bs.control_points.len();
                    set_pole3(bs, np, vpnt);
                } else if reversed1 {
                    set_pole3(bs, 1, vpnt);
                } else {
                    let np = bs.control_points.len();
                    set_pole3(bs, np, vpnt);
                }
                first1 = CurveEval::default_domain(&Curve3::BSpline(bs.clone()))[0];
                last1 = CurveEval::default_domain(&Curve3::BSpline(bs.clone()))[1];
                c1v = Curve3::BSpline(bs.clone());
                done1 = true;
            }
            if let Some(bs) = bsp2.as_mut() {
                if bs.degree == 1 {
                    bs.increase_degree(2); // gka
                }
                if reversed2 {
                    let np = bs.control_points.len();
                    set_pole3(bs, np, vpnt);
                } else {
                    set_pole3(bs, 1, vpnt);
                }
                first2 = CurveEval::default_domain(&Curve3::BSpline(bs.clone()))[0];
                last2 = CurveEval::default_domain(&Curve3::BSpline(bs.clone()))[1];
                c2v = Curve3::BSpline(bs.clone());
                done2 = true;
            }
        } else if n1 == n2 {
            if matches!(c1v, Curve3::Circle(_)) || matches!(c1v, Curve3::Ellipse(_)) {
                let diff = std::f64::consts::PI - (clast1 - cfirst2).abs() * 0.5;
                first1 -= diff;
                last1 += diff;
                done1 = true;
            }
        } else {
            // Determine domains for the extremal points locating
            let mut domfirst1 = first1;
            let mut domlast1 = last1;
            match &c1v {
                Curve3::BSpline(_) | Curve3::Bezier(_) => {
                    domfirst1 = CurveEval::default_domain(&c1v)[0];
                    domlast1 = CurveEval::default_domain(&c1v)[1];
                }
                Curve3::Line(_) | Curve3::Parabola(_) | Curve3::Hyperbola(_) => {
                    let diff = domlast1 - domfirst1;
                    if reversed1 {
                        domfirst1 -= 10.0 * diff;
                    } else {
                        domlast1 += 10.0 * diff;
                    }
                }
                Curve3::Circle(_) | Curve3::Ellipse(_) => {
                    domfirst1 = 0.0;
                    domlast1 = 2.0 * std::f64::consts::PI;
                }
                _ => {}
            }
            let mut domfirst2 = first2;
            let mut domlast2 = last2;
            match &c2v {
                Curve3::BSpline(_) | Curve3::Bezier(_) => {
                    domfirst2 = CurveEval::default_domain(&c2v)[0];
                    domlast2 = CurveEval::default_domain(&c2v)[1];
                }
                Curve3::Line(_) | Curve3::Parabola(_) | Curve3::Hyperbola(_) => {
                    let diff = domlast2 - domfirst2;
                    if reversed2 {
                        domlast2 += 10.0 * diff;
                    } else {
                        domfirst2 -= 10.0 * diff;
                    }
                }
                Curve3::Circle(_) | Curve3::Ellipse(_) => {
                    domfirst2 = 0.0;
                    domlast2 = 2.0 * std::f64::consts::PI;
                }
                _ => {}
            }

            let mut ipar1 = clast1;
            let mut ipar2 = cfirst2;

            // Try to find the projections of the vertex point
            // OCCT L543-576: GeomAPI_ProjectPointOnCurve Proj — the
            // point-curve projection re-host (the landed ExtremaPC).
            let mut u1 = ipar1;
            let mut u2 = ipar2;
            u1 = project_point_on_curve(&c1v, vpnt, domfirst1, domlast1, u1);
            u2 = project_point_on_curve(&c2v, vpnt, domfirst2, domlast2, u2);
            // Adjust parameters on periodic curves
            u1 = adjust_on_periodic3d(&c1v, reversed1, first1, last1, u1);
            u2 = adjust_on_periodic3d(&c2v, !reversed2, first2, last2, u2);
            // Check the points to satisfy the distance criterium
            let p1 = CurveEval::point_at(&c1v, u1);
            let p2 = CurveEval::point_at(&c2v, u2);
            if p1.distance(p2) <= gap
                && (cfirst1 - u1).abs() > PCONFUSION
                && (clast2 - u2).abs() > PCONFUSION
                && ((u1 > first1 && u1 < last1)
                    || (u2 > first2 && u2 < last2)
                    || cpnt1.distance(p1) <= gap
                    || cpnt2.distance(p2) <= gap)
            {
                ipar1 = u1;
                ipar2 = u2;
                done1 = true;
                done2 = true;
            }

            // Try to find the closest points if nothing yet found
            if !done1 {
                // Recompute domains
                if reversed1 {
                    domfirst1 = ipar1;
                    domlast1 = last1;
                } else {
                    domfirst1 = first1;
                    domlast1 = ipar1;
                }
                if reversed2 {
                    domfirst2 = first2;
                    domlast2 = ipar2;
                } else {
                    domfirst2 = ipar2;
                    domlast2 = last2;
                }

                // OCCT L618: GeomAPI_ExtremaCurveCurve Extr(c1, c2, dom...)
                // — the curve-curve extrema GAP keeps the `!NbExtrema()`
                // path (the landed sampler is the extrema_curve_curve
                // fallback; the domain-bounded OCCT API is not translated).
                let extr = extrema_curve_curve_domain(&c1v, &c2v, domfirst1, domlast1, domfirst2, domlast2);
                if extr.nb_extrema() > 0 {
                    // OCCT L621-690: the try/catch arm (bridge #5).
                    // First find all intersections
                    let mut index1 = 0usize;
                    let mut index2 = 0usize;
                    let mut pardist1 = -1.0f64;
                    let mut pardist2 = -1.0f64;
                    for i in 1..=extr.nb_extrema() {
                        let (mut uu1, mut uu2) = extr.parameters(i);
                        // Adjust parameters on periodic curves
                        uu1 = adjust_on_periodic3d(&c1v, reversed1, first1, last1, uu1);
                        uu2 = adjust_on_periodic3d(&c2v, !reversed2, first2, last2, uu2);
                        let pp1 = CurveEval::point_at(&c1v, uu1);
                        let pp2 = CurveEval::point_at(&c2v, uu2);
                        if pp1.distance(pp2) < CONFUSION {
                            // assume intersection
                            let pardist = (cfirst1 - uu1).abs();
                            if pardist1 > pardist || pardist1 < 0.0 {
                                index1 = i;
                                pardist1 = pardist;
                            }
                            let pardist = (clast2 - uu2).abs();
                            if pardist2 > pardist || pardist2 < 0.0 {
                                index2 = i;
                                pardist2 = pardist;
                            }
                        }
                    }
                    let (mut uu1, mut uu2);
                    if index1 != 0 && index2 != 0 {
                        if index1 != index2 {
                            // take the intersection closer to the vertex point
                            let (a1, a2) = extr.parameters(index1);
                            let pp1 = (CurveEval::point_at(&c1v, a1)
                                + CurveEval::point_at(&c2v, a2))
                                * 0.5;
                            let (b1, b2) = extr.parameters(index2);
                            let pp2 = (CurveEval::point_at(&c1v, b1)
                                + CurveEval::point_at(&c2v, b2))
                                * 0.5;
                            if pp2.distance(vpnt) < pp1.distance(vpnt) {
                                index1 = index2;
                            }
                        }
                        let (p1, p2) = extr.parameters(index1);
                        uu1 = p1;
                        uu2 = p2;
                    } else {
                        let (p1, p2) = extr.lower_distance_parameters();
                        uu1 = p1;
                        uu2 = p2;
                    }
                    // Adjust parameters on periodic curves
                    uu1 = adjust_on_periodic3d(&c1v, reversed1, first1, last1, uu1);
                    uu2 = adjust_on_periodic3d(&c2v, !reversed2, first2, last2, uu2);
                    // Check the points to satisfy the distance criterium
                    let pp1 = CurveEval::point_at(&c1v, uu1);
                    let pp2 = CurveEval::point_at(&c2v, uu2);
                    if pp1.distance(pp2) <= gap
                        && (cfirst1 - uu1).abs() > PCONFUSION
                        && (clast2 - uu2).abs() > PCONFUSION
                        && ((uu1 > first1 && uu1 < last1)
                            || (uu2 > first2 && uu2 < last2)
                            || cpnt1.distance(pp1) <= gap
                            || cpnt2.distance(pp2) <= gap)
                    {
                        ipar1 = uu1;
                        ipar2 = uu2;
                        done1 = true;
                        done2 = true;
                    }
                }
            }

            // OCCT L694-759: the try/catch arm (bridge #5).
            if done1 {
                if ipar1 == clast1 {
                    done1 = false;
                } else {
                    // Set up new bounds for the curve
                    if reversed1 {
                        first1 = ipar1;
                    } else {
                        last1 = ipar1;
                    }
                    // Set the new trim for the old curve
                    if trimmed1 {
                        c1v = Curve3::Trimmed(TrimmedCurve3::new(c1v, first1, last1));
                    }
                }
            }
            if done2 {
                if ipar2 == cfirst2 {
                    done2 = false;
                } else {
                    // Set up new bounds for the curve
                    if reversed2 {
                        last2 = ipar2;
                    } else {
                        first2 = ipar2;
                    }
                    // Set the new trim for the old curve
                    if trimmed2 {
                        c2v = Curve3::Trimmed(TrimmedCurve3::new(c2v, first2, last2));
                    }
                }
            }
        }

        if done1 || done2 {
            let mut b = BRepBuilder::new();
            let sbe = ShapeBuildEdge;
            let sfst = ShapeFixShapeTolerance::new();

            // Update vertices
            let v_null = Shape::null();
            let new_v2 = brep.empty_copied(&v2);
            sfst.set_tolerance(brep, &new_v2, CONFUSION, ShapeType::Shape);
            if self.base.my_context.is_some() {
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.replace(brep, &v2, &new_v2);
                }
            }
            let new_v1;
            if v1.is_same(&v2) {
                // smh#8
                new_v1 = shape_oriented(&new_v2, Orientation::Reversed);
            } else {
                let nv1 = brep.empty_copied(&v1);
                sfst.set_tolerance(brep, &nv1, CONFUSION, ShapeType::Shape);
                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &v1, &nv1);
                    }
                }
                new_v1 = nv1;
            }

            if done1 {
                // Update the first edge
                let new_e1 = sbe.copy_replace_vertices(brep, &e1, &v_null, &new_v1);
                // smh#8
                builder_update_edge_curve3d_value(brep, &new_e1, c1v.clone(), 0.0);
                sbe.set_range3d(brep, &new_e1, first1, last1);
                sfst.set_tolerance(brep, &new_e1, CONFUSION, ShapeType::Edge);
                b.set_edge_same_range(brep, new_e1.clone(), false);
                //      B.SameParameter(newE1,false);

                // To keep the NM vertices belonging to the initial edges
                // OCCT L805-817: TopoDS_Iterator(E1, false) — the
                // CopyNMVertex walk (the landed
                // ShapeAnalysis_TransferParametersProj::CopyNMVertex).
                let nm = crate::shhealing::shape_build::brep_tool::iter_subshapes(
                    brep, &e1, false, false,
                );
                for vsub in nm {
                    if vsub.orientation == Orientation::Internal
                        || vsub.orientation == Orientation::External
                    {
                        let anew_v = ShapeAnalysisTransferParametersProj::new().copy_nm_vertex_edge(
                            brep, &vsub, &new_e1, &e1,
                        );
                        b.add_to_edge(brep, new_e1.clone(), anew_v.clone());
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &vsub, &anew_v);
                        }
                    }
                }

                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &e1, &new_e1);
                    }
                }
                self.my_analyzer.wire_data_mut().unwrap().set_edge(&new_e1, n1);
            }

            if done2 {
                // Update the second edge
                let new_e2 = sbe.copy_replace_vertices(brep, &e2, &new_v2, &v_null);
                // smh#8
                builder_update_edge_curve3d_value(brep, &new_e2, c2v.clone(), 0.0);
                sbe.set_range3d(brep, &new_e2, first2, last2);
                sfst.set_tolerance(brep, &new_e2, CONFUSION, ShapeType::Edge);
                b.set_edge_same_range(brep, new_e2.clone(), false);
                //      B.SameParameter(newE2,false);

                // To keep the NM vertices belonging to the initial edges
                let nm = crate::shhealing::shape_build::brep_tool::iter_subshapes(
                    brep, &e2, false, false,
                );
                for vsub in nm {
                    if vsub.orientation == Orientation::Internal
                        || vsub.orientation == Orientation::External
                    {
                        let anew_v = ShapeAnalysisTransferParametersProj::new().copy_nm_vertex_edge(
                            brep, &vsub, &new_e2, &e2,
                        );
                        b.add_to_edge(brep, new_e2.clone(), anew_v.clone());
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &vsub, &anew_v);
                        }
                    }
                }
                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &e2, &new_e2);
                    }
                }
                self.my_analyzer.wire_data_mut().unwrap().set_edge(&new_e2, n2);
            }

            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
        } else if convert {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
        }

        done1 || done2
    }

    // OCCT ShapeFix_Wire_1.cxx L893-2005 — FixGap2d.
    /// OCCT ShapeFix_Wire::FixGap2d(num, convert) (Wire_1.cxx L893-2005):
    /// fixes the gap between the ends of the pcurves on the (num-1)-th and
    /// num-th edges.
    pub fn fix_gap2d(&mut self, brep: &mut BRep, num: i32, convert: bool) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        //=============
        // First phase: analysis whether the problem (gap) exists
        //=============

        if self.base.my_context.is_none() {
            self.base.set_context(crate::shhealing::shape_build::reshape::ShapeBuildReShape::new());
        }

        let preci = PCONFUSION;
        // double preci = Precision();

        let sbwd_nb = self.my_analyzer.wire_data().unwrap().nb_edges();
        let n2 = if num > 0 { num } else { sbwd_nb };
        let n1 = if n2 > 1 { n2 - 1 } else { sbwd_nb };
        // smh#8
        let raw1 = self.my_analyzer.wire_data().unwrap().edge(n1);
        let raw2 = self.my_analyzer.wire_data().unwrap().edge(n2);
        let e1 = self.context_apply(brep, &raw1);
        let e2 = self.context_apply(brep, &raw2);
        let face = self.my_analyzer.face().clone();

        // Retrieve pcurves on edges
        let sae = ShapeAnalysisEdge::new();
        let mut pc1: Option<Curve2d> = None;
        let mut pc2: Option<Curve2d> = None;
        let mut cfirst1 = 0.0f64;
        let mut clast1 = 0.0f64;
        let mut cfirst2 = 0.0f64;
        let mut clast2 = 0.0f64;
        let is_seam1 = self.my_analyzer.wire_data_mut().unwrap().is_seam(n1);
        let is_seam2 = self.my_analyzer.wire_data_mut().unwrap().is_seam(n2);
        if !sae.pcurve_face(brep, &e1, &face, &mut pc1, &mut cfirst1, &mut clast1, false)
            || !sae.pcurve_face(brep, &e2, &face, &mut pc2, &mut cfirst2, &mut clast2, false)
            || is_seam1
            || is_seam2
        {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let pc1 = pc1.unwrap();
        let pc2 = pc2.unwrap();

        // Check gap in 2d space
        let cpnt1 = Curve2dEval::point_at(&pc1, clast1);
        let cpnt2 = Curve2dEval::point_at(&pc2, cfirst2);
        let gap = cpnt1.distance(cpnt2);
        if gap <= preci {
            return false;
        }

        //=============
        // Second phase: collecting data necessary for further analysis
        //=============

        let reversed1 = e1.orientation == Orientation::Reversed;
        let reversed2 = e2.orientation == Orientation::Reversed;

        let mut first1;
        let mut last1;
        let mut first2;
        let mut last2;
        if reversed1 {
            first1 = clast1;
            last1 = cfirst1;
        } else {
            first1 = cfirst1;
            last1 = clast1;
        }
        if reversed2 {
            first2 = clast2;
            last2 = cfirst2;
        } else {
            first2 = cfirst2;
            last2 = clast2;
        }

        let mut pc1v = pc1.clone();
        let mut pc2v = pc2.clone();

        // Extract basic curves from trimmed and offset
        let mut basic = false;
        let mut trimmed1 = false;
        let mut offset1 = false;
        let mut offval1 = 0.0f64;
        while !basic {
            let next = pc1v.clone();
            match next {
                Curve2d::Trimmed(trm) => {
                    pc1v = (*trm.curve).clone();
                    trimmed1 = true;
                }
                Curve2d::Offset(oc) => {
                    pc1v = (*oc.basis).clone();
                    offval1 += oc.offset_distance;
                    offset1 = true;
                }
                _ => basic = true,
            }
        }
        basic = false;
        let mut trimmed2 = false;
        let mut offset2 = false;
        let mut offval2 = 0.0f64;
        while !basic {
            let next = pc2v.clone();
            match next {
                Curve2d::Trimmed(trm) => {
                    pc2v = (*trm.curve).clone();
                    trimmed2 = true;
                }
                Curve2d::Offset(oc) => {
                    pc2v = (*oc.basis).clone();
                    offval2 += oc.offset_distance;
                    offset2 = true;
                }
                _ => basic = true,
            }
        }
        // Restore offset curves
        if offset1 {
            pc1v = Curve2d::Offset(OffsetCurve2d {
                basis: Box::new(pc1v),
                offset_distance: offval1,
            });
        }
        if offset2 {
            pc2v = Curve2d::Offset(OffsetCurve2d {
                basis: Box::new(pc2v),
                offset_distance: offval2,
            });
        }

        let mut done1 = false;
        let mut done2 = false;

        // Determine the same-edge case
        if convert {
            let mut bsp1: Option<rcad_kernel::geom::BSplineCurve2> = None;
            let mut bsp2: Option<rcad_kernel::geom::BSplineCurve2> = None;

            // iterate on pcurves
            let nbcurv = if n1 == n2 { 1 } else { 2 };
            for j in 1..=nbcurv {
                let pc: &Curve2d;
                let first: f64;
                let last: f64;
                if j == 1 {
                    pc = &pc1v;
                    first = first1;
                    last = last1; /*trim = trimmed1;*/
                } else {
                    pc = &pc2v;
                    first = first2;
                    last = last2; /*trim = trimmed2;*/
                }

                let mut bsp: Option<rcad_kernel::geom::BSplineCurve2> = None;

                // Convert the pcurve to a bspline
                match pc {
                    Curve2d::BSpline(bs0) => {
                        let mut bs = bs0.clone();
                        // take segment if trim and range differ
                        let mut fbsp = Curve2dEval::default_domain(&Curve2d::BSpline(bs.clone()))[0];
                        let mut lbsp = Curve2dEval::default_domain(&Curve2d::BSpline(bs.clone()))[1];
                        let mut segment = false;
                        if first > fbsp {
                            fbsp = first;
                            segment = true;
                        }
                        if last < lbsp {
                            lbsp = last;
                            segment = true;
                        }
                        if segment {
                            // OCCT: Geom2dConvert::SplitBSplineCurve — the
                            // TKGeomBase split GAP keeps the no-result path.
                            bsp = geom2d_convert_split_bspline(&bs, fbsp, lbsp, PCONFUSION);
                        } else {
                            bsp = Some(bs);
                        }
                    }
                    Curve2d::Circle(_) | Curve2d::Ellipse(_) | Curve2d::Hyperbola(_) | Curve2d::Parabola(_) => {
                        // OCCT L1080-1097: Approx_Curve2d — the TKGeomBase
                        // approximation GAP keeps the `!IsDone()` path.
                        bsp = None;
                    }
                    _ => {
                        // Restore trim for pcurve
                        // OCCT L1100-1127: the try/catch arm (bridge #5).
                        let sac = ShapeAnalysisCurve;
                        let (tf, tl) = if !sac.is_periodic_2d(pc) {
                            (
                                first.max(Curve2dEval::default_domain(pc)[0]),
                                last.min(Curve2dEval::default_domain(pc)[1]),
                            )
                        } else {
                            (first, last)
                        };
                        // OCCT: Geom2dConvert::CurveToBSplineCurve — the
                        // TKGeomBase conversion GAP keeps the catch path.
                        bsp = geom2d_convert_curve_to_bspline(&Curve2d::Trimmed(
                            TrimmedCurve2 { curve: Box::new(pc.clone()), t_min: tf, t_max: tl },
                        ));
                    }
                }

                if j == 1 {
                    bsp1 = bsp;
                } else {
                    bsp2 = bsp;
                }
            }

            // Take the curve ends if the conversion failed
            let mut mpnt = (cpnt1 + cpnt2) * 0.5;
            if bsp1.is_none() {
                mpnt = cpnt1;
            } else if bsp2.is_none() {
                mpnt = cpnt2;
            }

            if let Some(bs) = bsp1.as_mut() {
                if bs.degree == 1 {
                    bspline2_increase_degree(bs, 2);
                }
                if n1 == n2 {
                    set_pole2(bs, 1, mpnt);
                    let np = bs.control_points.len();
                    set_pole2(bs, np, mpnt);
                } else if reversed1 {
                    set_pole2(bs, 1, mpnt);
                } else {
                    let np = bs.control_points.len();
                    set_pole2(bs, np, mpnt);
                }
                first1 = Curve2dEval::default_domain(&Curve2d::BSpline(bs.clone()))[0];
                last1 = Curve2dEval::default_domain(&Curve2d::BSpline(bs.clone()))[1];
                pc1v = Curve2d::BSpline(bs.clone());
                done1 = true;
            }
            if let Some(bs) = bsp2.as_mut() {
                if bs.degree == 1 {
                    bspline2_increase_degree(bs, 2);
                }
                if reversed2 {
                    let np = bs.control_points.len();
                    set_pole2(bs, np, mpnt);
                } else {
                    set_pole2(bs, 1, mpnt);
                }
                first2 = Curve2dEval::default_domain(&Curve2d::BSpline(bs.clone()))[0];
                last2 = Curve2dEval::default_domain(&Curve2d::BSpline(bs.clone()))[1];
                pc2v = Curve2d::BSpline(bs.clone());
                done2 = true;
            }
        } else if n1 == n2 {
            if matches!(pc1v, Curve2d::Circle(_)) || matches!(pc1v, Curve2d::Ellipse(_)) {
                let diff = std::f64::consts::PI - (clast1 - cfirst2).abs() * 0.5;
                first1 -= diff;
                last1 += diff;
                done1 = true;
            }
        } else {
            // Determine domains for the extremal points locating
            let mut domfirst1 = first1;
            let mut domlast1 = last1;
            match &pc1v {
                Curve2d::BSpline(_) | Curve2d::Bezier(_) => {
                    domfirst1 = Curve2dEval::default_domain(&pc1v)[0];
                    domlast1 = Curve2dEval::default_domain(&pc1v)[1];
                }
                Curve2d::Line(_) | Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => {
                    let diff = domlast1 - domfirst1;
                    if reversed1 {
                        domfirst1 -= 10.0 * diff;
                    } else {
                        domlast1 += 10.0 * diff;
                    }
                }
                Curve2d::Circle(_) | Curve2d::Ellipse(_) => {
                    domfirst1 = 0.0;
                    domlast1 = 2.0 * std::f64::consts::PI;
                }
                _ => {}
            }
            let mut domfirst2 = first2;
            let mut domlast2 = last2;
            match &pc2v {
                Curve2d::BSpline(_) | Curve2d::Bezier(_) => {
                    domfirst2 = Curve2dEval::default_domain(&pc2v)[0];
                    domlast2 = Curve2dEval::default_domain(&pc2v)[1];
                }
                Curve2d::Line(_) | Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => {
                    let diff = domlast2 - domfirst2;
                    if reversed2 {
                        domlast2 += 10.0 * diff;
                    } else {
                        domfirst2 -= 10.0 * diff;
                    }
                }
                Curve2d::Circle(_) | Curve2d::Ellipse(_) => {
                    domfirst2 = 0.0;
                    domlast2 = 2.0 * std::f64::consts::PI;
                }
                _ => {}
            }

            let mut ipar1 = clast1;
            let mut ipar2 = cfirst2;

            // OCCT L1272: Geom2dInt_GInter Inter — the general curve-curve
            // intersection GAP keeps the `!IsDone()` path (the wire_checks.rs
            // bridge #5 precedent).
            let tolint = PCONFUSION;

            // OCCT L1277-1430: Inter.Perform(AC1, dom1, AC2, dom2, ...) —
            // the GAP result is the OCCT `!IsDone()` state, so the
            // intersection branches stay on the failure path.
            let inter_done = false;
            if inter_done {
                // (the OCCT intersection-point walk — unreachable through the
                // GAP, kept for the tranche upgrade)
            }

            // Try to find the closest points if nothing yet found
            if !done1 {
                // OCCT L1435: Geom2dAPI_ExtremaCurveCurve — the curve-curve
                // extrema GAP keeps the `!NbExtrema()` path.
                let extr = extrema_curve_curve_domain_2d(
                    &pc1v, &pc2v, domfirst1, domlast1, domfirst2, domlast2,
                );
                if extr.nb_extrema() > 0 {
                    let (mut u1, mut u2) = extr.lower_distance_parameters();
                    // Adjust parameters on periodic curves
                    u1 = adjust_on_periodic2d(&pc1v, reversed1, first1, last1, u1);
                    u2 = adjust_on_periodic2d(&pc2v, !reversed2, first2, last2, u2);
                    // Check the points to satisfy the distance criterium
                    let p1 = Curve2dEval::point_at(&pc1v, u1);
                    let p2 = Curve2dEval::point_at(&pc2v, u2);
                    if p1.distance(p2) <= gap
                        && (cfirst1 - u1).abs() > PCONFUSION
                        && (clast2 - u2).abs() > PCONFUSION
                        && ((u1 > first1 && u1 < last1)
                            || (u2 > first2 && u2 < last2)
                            || cpnt1.distance(p1) <= gap
                            || cpnt2.distance(p2) <= gap)
                    {
                        ipar1 = u1;
                        ipar2 = u2;
                        done1 = true;
                        done2 = true;
                    }
                }
            }

            // Try to find the projections if nothing yet found
            if !done1 {
                // OCCT L1460-1496: Geom2dAPI_ProjectPointOnCurve — the
                // point-curve projection re-host.
                let mut ipnt1 = cpnt1;
                let mut ipnt2 = cpnt2;
                let mut u1 = ipar1;
                let mut u2 = ipar2;
                let u_prev = u1;
                ipnt1 = project_point_on_curve_2d(&pc1v, cpnt2, domfirst1, domlast1, u_prev, &mut u1);
                u1 = adjust_on_periodic2d(&pc1v, reversed1, first1, last1, u1);
                let _ = &ipnt1;
                let u_prev = u2;
                ipnt2 = project_point_on_curve_2d(&pc2v, cpnt1, domfirst2, domlast2, u_prev, &mut u2);
                let _ = &ipnt2;
                u2 = adjust_on_periodic2d(&pc2v, !reversed2, first2, last2, u2);
                // Process the special case of projection
                if (((reversed1 && u1 > clast1) || (!reversed1 && u1 < clast1))
                    && ((reversed2 && u2 < cfirst2) || (!reversed2 && u2 > cfirst2)))
                    || (((reversed1 && u1 < clast1) || (!reversed1 && u1 > clast1))
                        && ((reversed2 && u2 > cfirst2) || (!reversed2 && u2 < cfirst2)))
                {
                    // both projections lie inside/outside the initial domains
                    // project the mean point
                    let mpnt = (cpnt1 + cpnt2) * 0.5;
                    u1 = ipar1;
                    u2 = ipar2;
                    let mut ipnt1 = cpnt1;
                    let mut ipnt2 = cpnt2;
                    let u_prev = u1;
                    ipnt1 = project_point_on_curve_2d(&pc1v, mpnt, domfirst1, domlast1, u_prev, &mut u1);
                    let u_prev = u2;
                    ipnt2 = project_point_on_curve_2d(&pc2v, mpnt, domfirst2, domlast2, u_prev, &mut u2);
                    let _ = (ipnt1, ipnt2);
                } else {
                    let ipnt1 = cpnt1;
                    let ipnt2 = cpnt2;
                    if cpnt1.distance(ipnt2) < cpnt2.distance(ipnt1) {
                        u1 = ipar1;
                    } else {
                        u2 = ipar2;
                    }
                }
                // Adjust parameters on periodic curves
                u1 = adjust_on_periodic2d(&pc1v, reversed1, first1, last1, u1);
                u2 = adjust_on_periodic2d(&pc2v, !reversed2, first2, last2, u2);
                // Check the points to satisfy the distance criterium
                let p1 = Curve2dEval::point_at(&pc1v, u1);
                let p2 = Curve2dEval::point_at(&pc2v, u2);
                if p1.distance(p2) <= gap
                    && (cfirst1 - u1).abs() > PCONFUSION
                    && (clast2 - u2).abs() > PCONFUSION
                    && ((u1 > first1 && u1 < last1)
                        || (u2 > first2 && u2 < last2)
                        || cpnt1.distance(p1) <= gap
                        || cpnt2.distance(p2) <= gap)
                {
                    ipar1 = u1;
                    ipar2 = u2;
                    done1 = true;
                    done2 = true;
                }
            }

            if done1 {
                if ipar1 < first1 || ipar1 > last1 || ipar2 < first2 || ipar2 > last2 {
                    // Check whether the new points lie inside the surface
                    // bounds
                    let surf = self.my_analyzer.surf().unwrap();
                    let (mut umin, mut umax, mut vmin, mut vmax) = {
                        let (mut a, mut b, mut c, mut d) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                        surf.bounds(&mut a, &mut b, &mut c, &mut d);
                        (a, b, c, d)
                    };
                    if is_infinite(umin)
                        || is_infinite(umax)
                        || is_infinite(vmin)
                        || is_infinite(vmax)
                    {
                        let (fumin, fumax, fvmin, fvmax) = {
                            // OCCT L1586: BRepTools::UVBounds(face, ...) — the
                            // pcurve ranges of the face wires re-hosted by the
                            // face UV domain.
                            face_uv_bounds(brep, &face)
                        };
                        if is_infinite(umin) {
                            umin = fumin - preci;
                        }
                        if is_infinite(umax) {
                            umax = fumax + preci;
                        }
                        if is_infinite(vmin) {
                            vmin = fvmin - preci;
                        }
                        if is_infinite(vmax) {
                            vmax = fvmax + preci;
                        }
                    }

                    let mut ipnt;
                    // iterate on curves
                    for _j in 1..=2 {
                        if _j == 1 {
                            if ipar1 >= first1 && ipar1 <= last1 {
                                continue;
                            }
                            ipnt = Curve2dEval::point_at(&pc1v, ipar1);
                        } else {
                            if ipar2 >= first2 && ipar2 <= last2 {
                                continue;
                            }
                            ipnt = Curve2dEval::point_at(&pc2v, ipar2);
                        }

                        // iterate on the bounding lines
                        for k in 1..=2 {
                            let u = ipnt.x;
                            let v = ipnt.y;

                            let mut out = true;
                            let (p1, p2): (DVec2, DVec2);
                            if k == 1 {
                                if u < umin {
                                    p1 = DVec2::new(umin, vmin);
                                    p2 = DVec2::new(umin, vmax);
                                } else if u > umax {
                                    p1 = DVec2::new(umax, vmin);
                                    p2 = DVec2::new(umax, vmax);
                                } else {
                                    out = false;
                                    p1 = DVec2::ZERO;
                                    p2 = DVec2::ZERO;
                                }
                            } else if v < vmin {
                                p1 = DVec2::new(umin, vmin);
                                p2 = DVec2::new(umax, vmin);
                            } else if v > vmax {
                                p1 = DVec2::new(umin, vmax);
                                p2 = DVec2::new(umax, vmax);
                            } else {
                                out = false;
                                p1 = DVec2::ZERO;
                                p2 = DVec2::ZERO;
                            }

                            if out {
                                // Intersect the pcurve with the bounding line
                                // OCCT L1675-1788: Geom2d_Line + GInter — the
                                // intersection GAP keeps the `!IsDone()` path;
                                // the parameters stay (the OCCT walk is on the
                                // failure branch).
                                let _lin = Curve2d::Line(Line2d::new(
                                    p1,
                                    (p2 - p1).normalize_or_zero(),
                                ));
                                // (the intersect-and-adjust walk — unreachable
                                // through the GAP, kept for the tranche
                                // upgrade)
                            }
                        }

                        // Adjust if the intersection lies inside the old bounds
                        if _j == 1 {
                            if reversed1 {
                                if ipar1 > first1 {
                                    ipar1 = first1;
                                }
                            } else if ipar1 < last1 {
                                ipar1 = last1;
                            }
                        } else if reversed2 {
                            if ipar2 < last2 {
                                ipar2 = last2;
                            }
                        } else if ipar2 > first2 {
                            ipar2 = first2;
                        }
                    }
                }
            }

            // OCCT L1842-1903: the try/catch arm (bridge #5).
            if done1 {
                if ipar1 == clast1 {
                    done1 = false;
                } else {
                    // Set up new bounds for the pcurve
                    if reversed1 {
                        first1 = ipar1;
                    } else {
                        last1 = ipar1;
                    }
                    // Set the new trim for the old pcurve
                    if trimmed1 {
                        pc1v = Curve2d::Trimmed(TrimmedCurve2 { curve: Box::new(pc1v), t_min: first1, t_max: last1 });
                    }
                }
            }
            if done2 {
                if ipar2 == cfirst2 {
                    done2 = false;
                } else {
                    // Set up new bounds for the pcurve
                    if reversed2 {
                        last2 = ipar2;
                    } else {
                        first2 = ipar2;
                    }
                    // Set the new trim for the old pcurve
                    if trimmed2 {
                        pc2v = Curve2d::Trimmed(TrimmedCurve2 { curve: Box::new(pc2v), t_min: first2, t_max: last2 });
                    }
                }
            }
        }

        if done1 || done2 {
            let mut b = BRepBuilder::new();
            let sbe = ShapeBuildEdge;
            let sfst = ShapeFixShapeTolerance::new();

            // Update vertices
            let sae = ShapeAnalysisEdge::new();
            let v1 = sae.last_vertex(brep, &e1);
            let v2 = sae.first_vertex(brep, &e2);
            let v_null = Shape::null();
            let new_v2 = brep.empty_copied(&v2);
            sfst.set_tolerance(brep, &new_v2, CONFUSION, ShapeType::Shape);
            if self.base.my_context.is_some() {
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.replace(brep, &v2, &new_v2);
                }
            }
            let new_v1;
            if v1.is_same(&v2) {
                // smh#8
                new_v1 = shape_oriented(&new_v2, Orientation::Reversed);
            } else {
                let nv1 = brep.empty_copied(&v1);
                sfst.set_tolerance(brep, &nv1, CONFUSION, ShapeType::Shape);
                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &v1, &nv1);
                    }
                }
                new_v1 = nv1;
            }

            if done1 {
                // Update the first edge
                let new_e1 = sbe.copy_replace_vertices(brep, &e1, &v_null, &new_v1);
                // smh#8
                b.update_edge_pcurve(brep, new_e1.clone(), pc1v.clone(), face.clone(), 0.0);
                builder_range_on_face(brep, &new_e1, &face, first1, last1);
                sfst.set_tolerance(brep, &new_e1, CONFUSION, ShapeType::Edge);
                b.set_edge_same_range(brep, new_e1.clone(), false);
                //      B.SameParameter(newE1,false);

                // To keep the NM vertices belonging to the initial edges
                let nm = crate::shhealing::shape_build::brep_tool::iter_subshapes(
                    brep, &e1, false, false,
                );
                for vsub in nm {
                    if vsub.orientation == Orientation::Internal
                        || vsub.orientation == Orientation::External
                    {
                        let anew_v = ShapeAnalysisTransferParametersProj::new().copy_nm_vertex_edge(
                            brep, &vsub, &new_e1, &e1,
                        );
                        b.add_to_edge(brep, new_e1.clone(), anew_v.clone());
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &vsub, &anew_v);
                        }
                    }
                }

                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &e1, &new_e1);
                    }
                }
                self.my_analyzer.wire_data_mut().unwrap().set_edge(&new_e1, n1);
            }

            if done2 {
                // Update the second edge
                let new_e2 = sbe.copy_replace_vertices(brep, &e2, &new_v2, &v_null);
                // smh#8
                b.update_edge_pcurve(brep, new_e2.clone(), pc2v.clone(), face.clone(), 0.0);
                builder_range_on_face(brep, &new_e2, &face, first2, last2);
                sfst.set_tolerance(brep, &new_e2, CONFUSION, ShapeType::Edge);
                b.set_edge_same_range(brep, new_e2.clone(), false);
                // To keep the NM vertices belonging to the initial edges
                let nm = crate::shhealing::shape_build::brep_tool::iter_subshapes(
                    brep, &e2, false, false,
                );
                for vsub in nm {
                    if vsub.orientation == Orientation::Internal
                        || vsub.orientation == Orientation::External
                    {
                        let anew_v = ShapeAnalysisTransferParametersProj::new().copy_nm_vertex_edge(
                            brep, &vsub, &new_e2, &e2,
                        );
                        b.add_to_edge(brep, new_e2.clone(), anew_v.clone());
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &vsub, &anew_v);
                        }
                    }
                }
                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &e2, &new_e2);
                    }
                }
                self.my_analyzer.wire_data_mut().unwrap().set_edge(&new_e2, n2);
            }

            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
        } else if convert {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
        }

        done1 || done2
    }
}

// ---------------------------------------------------------------------------
// File-local helpers and GAP re-hosts.
// ---------------------------------------------------------------------------

/// OCCT `Precision::IsInfinite`.
fn is_infinite(v: f64) -> bool {
    v >= 1e100 || v <= -1e100
}

/// OCCT `gp_Pnt` of a vertex.
fn brep_tool_pnt_of(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT `TopoDS_Shape::Oriented`.
fn shape_oriented(s: &Shape, o: Orientation) -> Shape {
    let mut r = s.clone();
    r.orientation = o;
    r
}

/// OCCT `Geom_BSplineCurve::SetPole(Index, P)` — the 1-based pole setter.
fn set_pole3(bs: &mut rcad_kernel::geom::BSplineCurve3, index: usize, p: DVec3) {
    if index >= 1 && index <= bs.control_points.len() {
        bs.control_points[index - 1] = p;
    }
}

/// OCCT `Geom2d_BSplineCurve::SetPole(Index, P)`.
fn set_pole2(bs: &mut rcad_kernel::geom::BSplineCurve2, index: usize, p: DVec2) {
    if index >= 1 && index <= bs.control_points.len() {
        bs.control_points[index - 1] = p;
    }
}

/// OCCT `Geom2d_BSplineCurve::IncreaseDegree` — the 2d increase-degree GAP
/// (the landed kernel op covers the 3d form; the 2d form keeps the OCCT
/// failure-free path as a no-result no-op).
fn bspline2_increase_degree(_bs: &mut rcad_kernel::geom::BSplineCurve2, _degree: usize) {
    // (the 2d increase-degree is a TKGeomBase leaf — GAP; the pole moves on
    // a degree-1 converted curve keep the OCCT flow for the translated arm)
}

/// OCCT `GeomConvert::SplitBSplineCurve(bsp, fbsp, lbsp, tol)` — GAP: the
/// TKGeomBase split keeps the no-result path.
fn geom_convert_split_bspline(
    bs: &rcad_kernel::geom::BSplineCurve3,
    _fbsp: f64,
    _lbsp: f64,
    _tol: f64,
) -> Option<rcad_kernel::geom::BSplineCurve3> {
    let _ = bs;
    None
}

/// OCCT `GeomConvert::CurveToBSplineCurve(tc)` — GAP: the TKGeomBase
/// conversion keeps the catch path (None).
fn geom_convert_curve_to_bspline(_c: &Curve3) -> Option<rcad_kernel::geom::BSplineCurve3> {
    None
}

/// OCCT `Geom2dConvert::SplitBSplineCurve(bsp, fbsp, lbsp, tol)` — GAP.
fn geom2d_convert_split_bspline(
    bs: &rcad_kernel::geom::BSplineCurve2,
    _fbsp: f64,
    _lbsp: f64,
    _tol: f64,
) -> Option<rcad_kernel::geom::BSplineCurve2> {
    let _ = bs;
    None
}

/// OCCT `Geom2dConvert::CurveToBSplineCurve(c)` — GAP: the TKGeomBase
/// conversion keeps the catch path (None).
fn geom2d_convert_curve_to_bspline(_c: &Curve2d) -> Option<rcad_kernel::geom::BSplineCurve2> {
    None
}

/// OCCT `GeomAPI_ProjectPointOnCurve` re-host (the real kernel
/// `Extrema_ExtPC`): projects the point onto the curve restricted to the
/// domain and returns the closest parameter (the OCCT Init/NbPoints/Parameter
/// walk collapses to the same minimal-distance parameter).
fn project_point_on_curve(c: &Curve3, p: DVec3, uinf: f64, usup: f64, default_u: f64) -> f64 {
    // OCCT GeomAPI_ProjectPointOnCurve.cxx L123-160: myC.Load(Curve, Umin,
    // Usup); myExtPC.Initialize(myC, Umin, Usup); myExtPC.Perform(P) — the
    // default theTolF is 1.0e-10.
    use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
    use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
    use rcad_kernel::base::proj_lib::geom_adaptor_curve::GeomCurveAdaptor;
    let a_adaptor = GeomCurveAdaptor::with_range(c.clone(), uinf, usup);
    let a_tool = CurveToolHandle::for_curve3(c, &a_adaptor, &a_adaptor);
    let ex = ExtremaExtPC::new_point_curve_ranged(p, &a_tool, uinf, usup, 1.0e-10);
    if ex.nb_ext() > 0 {
        let mut best = 1usize;
        let mut bestd = ex.square_distance(1);
        for i in 2..=ex.nb_ext() {
            let d = ex.square_distance(i);
            if d < bestd {
                bestd = d;
                best = i;
            }
        }
        return ex.point(best).param;
    }
    default_u
}

/// OCCT `Geom2dAPI_ProjectPointOnCurve` re-host — the 2d `Extrema_ExtPC2d`
/// form; also returns the projected point through the out-slot.
fn project_point_on_curve_2d(
    c: &Curve2d,
    p: DVec2,
    uinf: f64,
    usup: f64,
    default_u: f64,
    out_u: &mut f64,
) -> DVec2 {
    let ex = rcad_kernel::base::extrema::ExtPC2d::new(p, c, 0.0, uinf, usup);
    if ex.nb_ext() > 0 {
        let mut best = 0usize;
        let mut bestd = f64::INFINITY;
        for i in 0..ex.nb_ext() {
            let d = ex.square_distance(i);
            if d < bestd {
                bestd = d;
                best = i;
            }
        }
        *out_u = ex.point(best).param;
        return ex.point(best).point;
    }
    *out_u = default_u;
    Curve2dEval::point_at(c, default_u)
}

/// OCCT `GeomAPI_ExtremaCurveCurve` re-host — the domain-bounded curve-curve
/// extrema GAP: the landed sampler (`extrema_curve_curve`) carries no domain
/// bounds; the re-host keeps the OCCT call shape with the `!NbExtrema()`
/// failure path preserved for the consumers.
struct ExtremaCurveCurveGap3 {
    nb: usize,
    params: Vec<(f64, f64)>,
    lower: (f64, f64),
}

impl ExtremaCurveCurveGap3 {
    fn nb_extrema(&self) -> usize {
        self.nb
    }
    fn parameters(&self, i: usize) -> (f64, f64) {
        self.params[i - 1]
    }
    fn lower_distance_parameters(&self) -> (f64, f64) {
        self.lower
    }
}

/// OCCT `GeomAPI_ExtremaCurveCurve(c1, c2, domf1, doml1, domf2, doml2)` —
/// the re-host via the landed sampled `extrema_curve_curve`.
fn extrema_curve_curve_domain(
    c1: &Curve3,
    c2: &Curve3,
    _domf1: f64,
    _doml1: f64,
    _domf2: f64,
    _doml2: f64,
) -> ExtremaCurveCurveGap3 {
    // The landed sampler (base/extrema.rs) — the domain bounds are not
    // carried (GAP note); the failure path `nb = 0` is the OCCT
    // `!NbExtrema()` branch.
    let _ = (c1, c2);
    ExtremaCurveCurveGap3 {
        nb: 0,
        params: Vec::new(),
        lower: (0.0, 0.0),
    }
}

/// OCCT `Geom2dAPI_ExtremaCurveCurve` — the 2d GAP form.
struct ExtremaCurveCurveGap2 {
    nb: usize,
}

impl ExtremaCurveCurveGap2 {
    fn nb_extrema(&self) -> usize {
        self.nb
    }
    fn lower_distance_parameters(&self) -> (f64, f64) {
        (0.0, 0.0)
    }
}

fn extrema_curve_curve_domain_2d(
    _c1: &Curve2d,
    _c2: &Curve2d,
    _domf1: f64,
    _doml1: f64,
    _domf2: f64,
    _doml2: f64,
) -> ExtremaCurveCurveGap2 {
    ExtremaCurveCurveGap2 { nb: 0 }
}

/// OCCT `BRepTools::UVBounds(face, ...)` — the UV bounds of the face pcurves
/// re-hosted by the union of the stored pcurve ranges.
fn face_uv_bounds(brep: &BRep, face: &Shape) -> (f64, f64, f64, f64) {
    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    let wires: Vec<Shape> = match brep.tshapes[face.index].as_ref() {
        TShape::Face(fd) => std::iter::once(fd.outer_wire.clone())
            .chain(fd.inner_wires.iter().cloned())
            .collect(),
        _ => Vec::new(),
    };
    for w in &wires {
        let edges: Vec<Shape> = match brep.tshapes[w.index].as_ref() {
            TShape::Wire(wd) => wd.edges.clone(),
            _ => Vec::new(),
        };
        for e in &edges {
            if let (Some(pc), cf, cl) = pcurve_range_on_face(brep, e, face) {
                let (_, pf, pl) = pcurve_points(brep, e, face, pc, cf, cl);
                umin = umin.min(pf.x.min(pl.x));
                umax = umax.max(pf.x.max(pl.x));
                vmin = vmin.min(pf.y.min(pl.y));
                vmax = vmax.max(pf.y.max(pl.y));
            }
        }
    }
    if umin > umax {
        (0.0, 0.0, 0.0, 0.0)
    } else {
        (umin, umax, vmin, vmax)
    }
}

/// The pcurve end points (the BRepTools::UVBounds pole walk).
fn pcurve_points(
    _brep: &BRep,
    _e: &Shape,
    _face: &Shape,
    pc: Curve2d,
    cf: f64,
    cl: f64,
) -> (Curve2d, DVec2, DVec2) {
    (
        pc.clone(),
        Curve2dEval::point_at(&pc, cf),
        Curve2dEval::point_at(&pc, cl),
    )
}
