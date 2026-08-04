//! IntPatch_ImpImpIntersection — intersection of two analytic surfaces.
//!
//! OCCT IntPatch_ImpImpIntersection.hxx / .cxx
//!
//! Handles all 15 pair combinations of Plane, Cylinder, Sphere, Cone, Torus
//! by converting to IntSurf_Quadric and dispatching to IntAna_QuadQuadGeo.

use super::GeomAbsSurfaceType;
use super::{
    AnaResultType, IntPatchIType, IntPatchLine, IntPatchPoint, IntPatchVertex, QuadQuadGeo, WLineType,
};
use crate::geomalgo::int_surf::quadric::Quadric;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3};

/// OCCT IntPatch_ImpImpIntersection.hxx L35-49: IntStatus
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntStatus {
    OK,
    InfiniteSectionCurve,
    Fail,
}

/// IntPatch_ImpImpIntersection
///
/// Fields (L117-124):
///   myDone: IntStatus
///   empt, tgte, oppo: bool
///   spnt: Vec<IntPatchPoint>
///   slin: Vec<IntPatchLine>
///   solrst: IntPatch_TheSOnBounds (skipped — rcad does boundary clipping elsewhere)
pub struct ImpImpIntersection {
    my_done: IntStatus,
    empt: bool,
    tgte: bool,
    oppo: bool,
    spnt: Vec<IntPatchPoint>,
    slin: Vec<IntPatchLine>,
    /// OCCT L2536: SameSurf — the two surfaces are coincident.
    same_surf: bool,
}

impl ImpImpIntersection {
    pub fn new() -> Self {
        Self {
            my_done: IntStatus::Fail,
            empt: true,
            tgte: false,
            oppo: false,
            spnt: Vec::new(),
            slin: Vec::new(),
            same_surf: false,
        }
    }

    pub fn is_done(&self) -> bool {
        self.my_done == IntStatus::OK
    }
    pub fn status(&self) -> IntStatus {
        self.my_done
    }
    pub fn is_empty(&self) -> bool {
        self.empt
    }
    pub fn tangent_faces(&self) -> bool {
        self.tgte
    }
    pub fn opposite_faces(&self) -> bool {
        self.oppo
    }
    pub fn nb_lines(&self) -> usize {
        self.slin.len()
    }
    pub fn line(&self, i: usize) -> &IntPatchLine {
        &self.slin[i]
    }
    pub fn line_mut(&mut self, i: usize) -> &mut IntPatchLine {
        &mut self.slin[i]
    }
    pub fn slin_ref(&self) -> &[IntPatchLine] {
        &self.slin
    }
    pub fn nb_points(&self) -> usize {
        self.spnt.len()
    }
    pub fn point(&self, i: usize) -> &IntPatchPoint {
        &self.spnt[i]
    }

    // =====================================================================
    // OCCT L73-79: Perform(S1, D1, S2, D2, TolArc, TolTang, theIsReqToKeepRLine)
    // rcad: Surface3 instead of Adaptor3d_Surface; the UV domains D1/D2 (the
    // corrected FF UV rectangles) replace the TopolTool.
    // =====================================================================
    pub fn perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        uv1: [f64; 4],
        uv2: [f64; 4],
        tol_arc: f64,
        tol_tang: f64,
    ) {
        // OCCT L2525-2533: myDone=Fail, spnt/slin clear, empt/tgte/oppo init
        self.my_done = IntStatus::Fail;
        self.spnt.clear();
        self.slin.clear();
        self.empt = true;
        self.tgte = false;
        self.oppo = false;
        self.same_surf = false;

        // OCCT L2529: isPostProcessingRequired = true
        let mut is_post_processing_required = true;
        // OCCT L2538: multpoint = false
        let mut multpoint = false;

        // OCCT L2548-2556: SetQuad — convert Surface3 to Quadric + type index
        let Some(q1) = Quadric::from_surface3(s1) else {
            return;
        };
        let Some(q2) = Quadric::from_surface3(s2) else {
            return;
        };
        // OCCT L2555: typs1, typs2 — surface type enums (rcad: inferred from Quadric)
        let typs1 = q1.surface_type();
        let typs2 = q2.surface_type();
        // OCCT L2553: bool bEmpty = false
        let mut b_empty = false;

        // OCCT L2558-2562: if (!iT1 || !iT2) throw ConstructionError
        let i_t1 = quad_type_index(&q1);
        let i_t2 = quad_type_index(&q2);
        if i_t1 == 0 || i_t2 == 0 {
            return;
        }

        // OCCT L2564-2565: bReverse = iT1 > iT2, iTT = iT1 * 10 + iT2
        let b_reverse = i_t1 > i_t2;
        let i_tt = i_t1 * 10 + i_t2;

        // OCCT L2567-2769: switch(iTT)
        match i_tt {
            // OCCT L2569-2574: case 11 Plane/Plane: if (!IntPP(...)) return;
            11 => {
                if !self.int_pp(&q1, &q2, tol_tang) {
                    return;
                }
            }
            // OCCT L2577-2594: case 12/21 Plane/Cylinder: H from surface bounds
            12 | 21 => {
                // OCCT L2579-2586: H = VMax - VMin of the cylinder surface V domain
                let a_s_cyl_uv = if b_reverse { uv1 } else { uv2 };
                let v_min = a_s_cyl_uv[2];
                let v_max = a_s_cyl_uv[3];
                let h = if rcad_kernel::precision::is_negative_infinite_value(v_min)
                    || rcad_kernel::precision::is_positive_infinite_value(v_max)
                {
                    0.0
                } else {
                    v_max - v_min
                };
                // OCCT L2588-2591: if (!IntPCy(...)) return;
                if !self.int_pcy(&q1, &q2, TOL_ANG, tol_tang, b_reverse, h) {
                    return;
                }
                // OCCT L2592: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2596-2604: case 13/31 Plane/Cone
            13 | 31 => {
                // OCCT L2598-2601: if (!IntPCo(...)) return;
                if !self.int_pco(&q1, &q2, TOL_ANG, tol_tang, b_reverse, &mut multpoint) {
                    return;
                }
                // OCCT L2602: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2606-2614: case 14/41 Plane/Sphere
            14 | 41 => {
                // OCCT L2608-2611: if (!IntPSp(...)) return;
                if !self.int_psp(&q1, &q2, TOL_ANG, tol_tang, b_reverse) {
                    return;
                }
                // OCCT L2612: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2616-2624: case 15/51 Plane/Torus
            15 | 51 => {
                // OCCT L2618-2621: if (!IntPTo(...)) return;
                if !self.int_pto(&q1, &q2, tol_tang, b_reverse) {
                    return;
                }
                // OCCT L2622: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2626-2677: case 22 Cylinder/Cylinder (aBox1,aBox2,a2DTol)
            22 => {
                // OCCT L2649-2652: a2DTol = min(1e-4, min(S1->UResolution(TolTang),
                // S2->UResolution(TolTang))).
                let a_2d_tol = 1.0e-4_f64
                    .min(rcad_kernel::topo::topods::u_resolution_for_surface(s1, tol_tang))
                    .min(rcad_kernel::topo::topods::u_resolution_for_surface(s2, tol_tang));
                // OCCT L2657-2658: myDone = IntCyCy(...); if Fail return.
                self.int_cycy(&q1, &q2, tol_tang, uv1, uv2, a_2d_tol, &mut multpoint);
                if self.my_done == IntStatus::Fail {
                    return;
                }
                // OCCT L2665: bEmpty = empt
                b_empty = self.empt;
                // OCCT L2666-2674: no geometric solution (numeric WLine path)
                // -> skip the post-processing.
                if !self.slin.is_empty() && self.slin[0].is_wline() {
                    is_post_processing_required = false;
                }
            }
            // OCCT L2679-2687: case 23/32 Cylinder/Cone
            23 | 32 => {
                // OCCT L2681-2684: if (!IntCyCo(...)) return;
                if !self.int_cyco(&q1, &q2, tol_tang, b_reverse, &mut multpoint) {
                    return;
                }
                // OCCT L2685: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2689-2697: case 24/42 Cylinder/Sphere
            24 | 42 => {
                // OCCT L2691-2694: if (!IntCySp(...)) return;
                if !self.int_cysp(&q1, &q2, tol_tang, b_reverse, &mut multpoint) {
                    return;
                }
                // OCCT L2695: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2699-2707: case 25/52 Cylinder/Torus
            25 | 52 => {
                // OCCT L2701-2704: if (!IntCyTo(...)) return;
                if !self.int_cyto(&q1, &q2, tol_tang, b_reverse) {
                    return;
                }
                // OCCT L2705: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2709-2716: case 33 Cone/Cone
            33 => {
                // OCCT L2710-2713: if (!IntCoCo(...)) return;
                if !self.int_coco(&q1, &q2, tol_tang, &mut multpoint) {
                    return;
                }
                // OCCT L2714: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2718-2726: case 34/43 Cone/Sphere
            34 | 43 => {
                // OCCT L2720-2723: if (!IntCoSp(...)) return;
                if !self.int_cosp(&q1, &q2, tol_tang, b_reverse, &mut multpoint) {
                    return;
                }
                // OCCT L2724: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2728-2735: case 35/53 Cone/Torus (no bEmpty)
            35 | 53 => {
                // OCCT L2730-2733: if (!IntCoTo(...)) return;
                if !self.int_coto(&q1, &q2, tol_tang, b_reverse) {
                    return;
                }
            }
            // OCCT L2737-2744: case 44 Sphere/Sphere
            44 => {
                // OCCT L2738-2741: if (!IntSpSp(...)) return;
                if !self.int_spsp(&q1, &q2, tol_tang) {
                    return;
                }
                // OCCT L2742: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2746-2754: case 45/54 Sphere/Torus
            45 | 54 => {
                // OCCT L2748-2751: if (!IntSpTo(...)) return;
                if !self.int_spto(&q1, &q2, tol_tang, b_reverse) {
                    return;
                }
                // OCCT L2752: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2756-2763: case 55 Torus/Torus
            55 => {
                // OCCT L2757-2760: if (!IntToTo(...)) return;
                if !self.int_toto(&q1, &q2, tol_tang) {
                    return;
                }
                // OCCT L2761: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2765-2768: default throw ConstructionError
            _ => {
                return;
            }
        }

        // OCCT L2771-2779: if bEmpty { myDone = OK; return }
        if b_empty {
            if self.my_done == IntStatus::Fail {
                self.my_done = IntStatus::OK;
            }
            return;
        }

        // OCCT L2782-2934: isPostProcessingRequired block
        if is_post_processing_required {
        let same_surf = self.same_surf;
        let mut all1 = false;
        let mut all2 = false;
        let mut nosolon_s1 = false;
        let mut nosolon_s2 = false;
        let mut edg1: Vec<super::so_on_bounds::Segment> = Vec::new();
        let mut edg2: Vec<super::so_on_bounds::Segment> = Vec::new();
        let mut pnt1: Vec<super::so_on_bounds::PathPoint> = Vec::new();
        let mut pnt2: Vec<super::so_on_bounds::PathPoint> = Vec::new();

        let mut a_func = super::so_on_bounds::ArcFunction::new();
        let mut solrst = super::so_on_bounds::SOnBounds::new();
        // OCCT D1/D2 (Adaptor3d_TopolTool) — created here for the UV-rectangle
        // domain; kept alive for the PutPointsOnLine call below.
        let mut d1: Option<super::so_on_bounds::Domain> = None;
        let mut d2: Option<super::so_on_bounds::Domain> = None;

        if !same_surf {
            // OCCT L2786-2787: AFunc.SetQuadric(quad2); AFunc.Set(S1);
            a_func.set_quadric(q2.clone());
            a_func.set_surface(s1.clone());
            // OCCT L2789: solrst.Perform(AFunc, D1, TolArc, TolTang);
            d1 = Some(super::so_on_bounds::Domain::new(uv1[0], uv1[1], uv1[2], uv1[3]));
            solrst.perform(&mut a_func, d1.as_mut().unwrap(), tol_arc, tol_tang, false);
            if !solrst.is_done() {
                return;
            }
            // OCCT L2795-2798: AllArcSolution && typs1 == typs2 -> all1 = true
            if solrst.all_arc_solution() && typs1 == typs2 {
                all1 = true;
            }
            let nbpt = solrst.nb_points();
            let nbseg = solrst.nb_segments();
            for i in 1..=nbpt {
                let a_pt = solrst.point(i).clone();
                pnt1.push(a_pt);
            }
            for i in 1..=nbseg {
                let a_segm = solrst.segment(i).clone();
                edg1.push(a_segm);
            }
            nosolon_s1 = (nbpt == 0) && (nbseg == 0);
            if nosolon_s1 && all1 {
                // case of a face without restrictions.
                all1 = false;
            }
        } else {
            nosolon_s1 = true;
        }

        if !same_surf {
            // OCCT L2825-2826: AFunc.SetQuadric(quad1); AFunc.Set(S2);
            a_func.set_quadric(q1.clone());
            a_func.set_surface(s2.clone());
            // OCCT L2828: solrst.Perform(AFunc, D2, TolArc, TolTang);
            d2 = Some(super::so_on_bounds::Domain::new(uv2[0], uv2[1], uv2[2], uv2[3]));
            solrst.perform(&mut a_func, d2.as_mut().unwrap(), tol_arc, tol_tang, false);
            if !solrst.is_done() {
                return;
            }
            if solrst.all_arc_solution() && typs1 == typs2 {
                all2 = true;
            }
            let nbpt = solrst.nb_points();
            let nbseg = solrst.nb_segments();
            for i in 1..=nbpt {
                let a_pt = solrst.point(i).clone();
                pnt2.push(a_pt);
            }
            for i in 1..=nbseg {
                let a_segm = solrst.segment(i).clone();
                edg2.push(a_segm);
            }
            nosolon_s2 = (nbpt == 0) && (nbseg == 0);
            if nosolon_s2 && all2 {
                // case of a face without restrictions.
                all2 = false;
            }
        } else {
            nosolon_s2 = true;
        }

        if same_surf || (all1 && all2) {
            // faces "paralleles" parfaites
            self.empt = false;
            self.tgte = true;
            self.slin.clear();
            self.spnt.clear();

            let ptreference;
            match typs1 {
                crate::geomalgo::int_surf::quadric::QuadricType::Plane => {
                    ptreference = q1.location();
                }
                crate::geomalgo::int_surf::quadric::QuadricType::Cylinder => {
                    ptreference = q1.value(0.0, 0.0);
                }
                crate::geomalgo::int_surf::quadric::QuadricType::Sphere => {
                    ptreference = q1.value(std::f64::consts::FRAC_PI_4, std::f64::consts::FRAC_PI_4);
                }
                crate::geomalgo::int_surf::quadric::QuadricType::Cone => {
                    ptreference = q1.value(0.0, 10.0);
                }
                crate::geomalgo::int_surf::quadric::QuadricType::Torus => {
                    ptreference = q1.value(0.0, 0.0);
                }
                _ => {
                    ptreference = DVec3::ZERO;
                }
            }

            self.oppo = q1.normale(ptreference).dot(q2.normale(ptreference)) < 0.0;
            self.my_done = IntStatus::OK;
            return;
        }

        if !nosolon_s1 || !nosolon_s2 {
            self.empt = false;
            // OCCT L2910: PutPointsOnLine(S1, S2, pnt1, slin, true, D1, quad1, quad2, multpoint, TolArc);
            super::restriction::put_points_on_line(
                s1, s2, &pnt1, &mut self.slin, true, d1.as_ref().unwrap(), &q1, &q2, multpoint,
                tol_arc,
            );
            // OCCT L2912: PutPointsOnLine(S1, S2, pnt2, slin, false, D2, quad2, quad1, multpoint, TolArc);
            super::restriction::put_points_on_line(
                s1, s2, &pnt2, &mut self.slin, false, d2.as_ref().unwrap(), &q2, &q1, multpoint,
                tol_arc,
            );

            if !edg1.is_empty() {
                super::restriction::process_segments(&edg1, &mut self.slin, &q1, &q2, true, tol_arc);
            }
            if !edg2.is_empty() {
                super::restriction::process_segments(&edg2, &mut self.slin, &q1, &q2, false, tol_arc);
            }
            if !edg1.is_empty() || !edg2.is_empty() {
                // OCCT L2927: ProcessRLine(slin, quad1, quad2, TolArc, theIsReqToKeepRLine);
                super::restriction::process_r_line(&mut self.slin, &q1, &q2, tol_arc, false);
            }
        } else {
            self.empt = self.slin.is_empty() && self.spnt.is_empty();
        }
        }

        // OCCT L2976-2995: ComputeVertexParameters(TolArc) for each GLine
        // (IntPatch_GLine.cxx L421-1090: filter, sort, dedup), for each ALine
        // (IntPatch_ALine.cxx L77-679: filter, sort, dedup) and for each RLine
        // (IntPatch_RLine.cxx L143-434: filter, sort, dedup).
        for line in self.slin.iter_mut() {
            let is_gline = matches!(
                line.line_type,
                IntPatchIType::Line
                    | IntPatchIType::Circle
                    | IntPatchIType::Ellipse
                    | IntPatchIType::Parabola
                    | IntPatchIType::Hyperbola
            );
            if is_gline {
                line.compute_vertex_parameters_gline();
            } else if line.line_type == IntPatchIType::Analytic {
                // OCCT L2985-2989: IntPatch_ALine::ComputeVertexParameters
                // (IntPatch_ALine.cxx L77-679) sorts the ALine vertices by
                // parameter so that IntPatch_ALineToWLine::MakeWLine walks them
                // in order.  The ALine's vertices live on the wrapped
                // IntAnaCurve.
                if let Some(ac) = line.a_curve.as_mut() {
                    ac.compute_vertex_parameters_aline(tol_arc);
                }
            } else if line.line_type == IntPatchIType::Restriction {
                // OCCT L2990-2994: IntPatch_RLine::ComputeVertexParameters.
                line.compute_vertex_parameters_rline();
            }
        }

        // OCCT L2997-3052: place 2 vertices on the GLine curves (Circle/
        // Ellipse) that have none (ElCLib::Value at param 0 and 2π; the 3D
        // point is the same for both — only the line parameter differs).
        for line in self.slin.iter_mut() {
            match line.line_type {
                IntPatchIType::Circle => {
                    if line.nb_vertex() == 0 {
                        if let Curve3::Circle(circ) = &line.curve {
                            let p = super::elclib::circle_value(circ, 0.0);
                            let (u1, v1) = q1.parameters(p);
                            let (u2, v2) = q2.parameters(p);
                            let mut point = IntPatchVertex::default();
                            point.set_value(p, tol_arc, false);
                            point.set_parameters(u1, v1, u2, v2);
                            point.set_parameter(0.0);
                            line.add_vertex(point.clone());
                            point.set_value(p, tol_arc, false);
                            point.set_parameters(u1, v1, u2, v2);
                            point.set_parameter(std::f64::consts::TAU);
                            line.add_vertex(point);
                        }
                    }
                }
                IntPatchIType::Ellipse => {
                    if line.nb_vertex() == 0 {
                        if let Curve3::Ellipse(e) = &line.curve {
                            let p = super::elclib::ellipse_value(e, 0.0);
                            let (u1, v1) = q1.parameters(p);
                            let (u2, v2) = q2.parameters(p);
                            let mut point = IntPatchVertex::default();
                            point.set_value(p, tol_arc, false);
                            point.set_parameters(u1, v1, u2, v2);
                            point.set_parameter(0.0);
                            line.add_vertex(point.clone());
                            point.set_value(p, tol_arc, false);
                            point.set_parameters(u1, v1, u2, v2);
                            point.set_parameter(std::f64::consts::TAU);
                            line.add_vertex(point);
                        }
                    }
                }
                _ => {}
            }
        }
    }

    // =====================================================================
    // 15 IntXX() functions — each delegates to the int_xx.rs 1:1 translation
    // =====================================================================

    fn int_pp(&mut self, q1: &Quadric, q2: &Quadric, tol_tang: f64) -> bool {
        // OCCT L2570: if (!IntPP(quad1, quad2, Tolang, TolTang, SameSurf, slin))
        // return;  IntPP has no Empty out-param (plane-plane) — empt keeps its
        // initial value and is refined by the SameSurf / SOnBounds blocks.
        let done = super::int_xx::int_pp(
            q1,
            q2,
            1e-8,
            tol_tang,
            &mut self.same_surf,
            &mut self.slin,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L3157-3345: IntPCy — Plane/Cylinder intersection.
    // Delegates to the 1:1 translation (int_xx.rs) which handles the tangent
    // line (NbSol==1), the two-line transitions, AdjustToSeam for the Circle,
    // and the Empty/OK flags exactly as OCCT IntPCy.
    fn int_pcy(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol_ang: f64,
        tol: f64,
        b_reverse: bool,
        h: f64,
    ) -> bool {
        // OCCT L2588: if (!IntPCy(quad1, quad2, Tolang, TolTang, bReverse,
        // empt, slin, H)) return;
        let done = super::int_xx::int_pcy(
            q1,
            q2,
            tol_ang,
            tol,
            b_reverse,
            &mut self.empt,
            &mut self.slin,
            h,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L3446-3706: IntPCo — Plane/Cone intersection.  Delegates to the 1:1
    // translation (int_xx.rs) which handles the apex tangent line, the two-line
    // transitions (Multpoint=true) and the Empty/OK flags exactly as OCCT.
    fn int_pco(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol_ang: f64,
        tol: f64,
        b_reverse: bool,
        multpoint: &mut bool,
    ) -> bool {
        // OCCT L2598: if (!IntPCo(...)) return;
        let done = super::int_xx::int_pco(
            q1,
            q2,
            tol_ang,
            tol,
            b_reverse,
            &mut self.empt,
            multpoint,
            &mut self.slin,
            &mut self.spnt,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L3352-3432: IntPSp — Plane/Sphere intersection.  Delegates to the
    // 1:1 translation (int_xx.rs).
    fn int_psp(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol_ang: f64,
        tol: f64,
        b_reverse: bool,
    ) -> bool {
        // OCCT L2608: if (!IntPSp(...)) return;
        let done = super::int_xx::int_psp(
            q1,
            q2,
            tol_ang,
            tol,
            b_reverse,
            &mut self.empt,
            &mut self.slin,
            &mut self.spnt,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L3708-3800: IntPTo — Plane/Torus intersection.  Delegates to the 1:1
    // translation (int_xx.rs).
    fn int_pto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) -> bool {
        // OCCT L2618: if (!IntPTo(...)) return;
        let done = super::int_xx::int_pto(
            q1,
            q2,
            tol,
            b_reverse,
            &mut self.empt,
            &mut self.slin,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    /// OCCT case 22 (L2626-2677): Cylinder/Cylinder.  Delegates to the 1:1
    /// IntCyCy translation (int_cycy.rs) which handles both the analytic
    /// (CyCyAnalyticalIntersect) and the numeric (CyCyNoGeometric) paths.
    fn int_cycy(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol: f64,
        uv1: [f64; 4],
        uv2: [f64; 4],
        tol_2d: f64,
        multpoint: &mut bool,
    ) {
        self.my_done = super::int_cycy::int_cycy(
            q1,
            q2,
            tol,
            tol_2d,
            uv1,
            uv2,
            &mut self.empt,
            &mut self.same_surf,
            multpoint,
            &mut self.slin,
            &mut self.spnt,
        );
    }

    // OCCT L8465-8610: IntCyCo — Cylinder/Cone intersection.  Delegates to the
    // 1:1 translation (int_xx.rs) which handles the NoGeometricSolution
    // IntAna_IntQuadQuad fallback internally.
    fn int_cyco(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol: f64,
        b_reverse: bool,
        multpoint: &mut bool,
    ) -> bool {
        // OCCT L2681: if (!IntCyCo(...)) return;
        let done = super::int_xx::int_cyco(
            q1,
            q2,
            tol,
            b_reverse,
            &mut self.empt,
            multpoint,
            &mut self.slin,
            &mut self.spnt,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L8263-8437: IntCySp — Cylinder/Sphere intersection.  Delegates to
    // the 1:1 translation (int_xx.rs).
    fn int_cysp(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol: f64,
        b_reverse: bool,
        multpoint: &mut bool,
    ) -> bool {
        // OCCT L2691: if (!IntCySp(...)) return;
        let done = super::int_xx::int_cysp(
            q1,
            q2,
            tol,
            b_reverse,
            &mut self.empt,
            multpoint,
            &mut self.slin,
            &mut self.spnt,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L8072-8240: IntCyTo — Cylinder/Torus intersection.  Delegates to the
    // 1:1 translation (int_xx.rs).
    fn int_cyto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) -> bool {
        // OCCT L2701: if (!IntCyTo(...)) return;
        let done = super::int_xx::int_cyto(
            q1,
            q2,
            tol,
            b_reverse,
            &mut self.empt,
            &mut self.slin,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L9022-9260: IntCoCo — Cone/Cone intersection.  Delegates to the 1:1
    // translation (int_xx.rs).
    fn int_coco(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol: f64,
        multpoint: &mut bool,
    ) -> bool {
        // OCCT L2710: if (!IntCoCo(...)) return;
        let done = super::int_xx::int_coco(
            q1,
            q2,
            tol,
            &mut self.empt,
            &mut self.same_surf,
            multpoint,
            &mut self.slin,
            &mut self.spnt,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT L9349-9590: IntCoSp — Cone/Sphere intersection.  Delegates to the
    // 1:1 translation (int_xx.rs).
    fn int_cosp(
        &mut self,
        q1: &Quadric,
        q2: &Quadric,
        tol: f64,
        b_reverse: bool,
        multpoint: &mut bool,
    ) -> bool {
        // OCCT L2720: if (!IntCoSp(...)) return;
        let done = super::int_xx::int_cosp(
            q1,
            q2,
            tol,
            b_reverse,
            &mut self.empt,
            multpoint,
            &mut self.slin,
            &mut self.spnt,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT IntCoTo — Cone/Torus intersection.  Delegates to the 1:1 translation
    // (int_xx.rs).
    fn int_coto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) -> bool {
        // OCCT L2730: if (!IntCoTo(...)) return;
        let done = super::int_xx::int_coto(q1, q2, tol, b_reverse, &mut self.empt, &mut self.slin);
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT IntSpSp — Sphere/Sphere intersection.  Delegates to the 1:1
    // translation (int_xx.rs).
    fn int_spsp(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) -> bool {
        // OCCT L2738: if (!IntSpSp(...)) return;
        let done = super::int_xx::int_spsp(
            q1,
            q2,
            tol,
            &mut self.empt,
            &mut self.same_surf,
            &mut self.slin,
            &mut self.spnt,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT IntSpTo — Sphere/Torus intersection.  Delegates to the 1:1
    // translation (int_xx.rs).
    fn int_spto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) -> bool {
        // OCCT L2748: if (!IntSpTo(...)) return;
        let done = super::int_xx::int_spto(q1, q2, tol, b_reverse, &mut self.empt, &mut self.slin);
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }

    // OCCT IntToTo — Torus/Torus intersection.  Delegates to the 1:1
    // translation (int_xx.rs).
    fn int_toto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) -> bool {
        // OCCT L2757: if (!IntToTo(...)) return;
        let done = super::int_xx::int_toto(
            q1,
            q2,
            tol,
            &mut self.same_surf,
            &mut self.empt,
            &mut self.slin,
        );
        if !done {
            return false;
        }
        self.my_done = IntStatus::OK;
        true
    }
}

/// OCCT: SetQuad equivalent — return type index (1=Plane, 2=Cyl, 3=Cone, 4=Sphere, 5=Torus)
fn quad_type_index(q: &Quadric) -> i32 {
    match super::geom_abs_of_quadric(q.surface_type()) {
        GeomAbsSurfaceType::Plane => 1,
        GeomAbsSurfaceType::Cylinder => 2,
        GeomAbsSurfaceType::Cone => 3,
        GeomAbsSurfaceType::Sphere => 4,
        GeomAbsSurfaceType::Torus => 5,
        _ => 0,
    }
}

/// Virtual tolerance angle used in Place/Cylinder dispatch
/// OCCT ProcessBounds (IntPatch_ImpImpIntersection.cxx L4683-4838).
///
/// Adds the endpoints (theFPar/theLPar) of the current ALine as IntPatch_Point
/// vertices.  When an endpoint coincides (within theTol) with a vertex of a
/// previously stored ALine, that vertex is reused (marked multiple) instead of
/// creating a new one.
pub(crate) fn process_bounds(
    alig: &mut super::int_quad_quad::IntAnaCurve,
    slin: &[IntPatchLine],
    quad1: &Quadric,
    quad2: &Quadric,
    procf: &mut bool,
    ptf: DVec3,
    first: f64,
    procl: &mut bool,
    ptl: DVec3,
    last: f64,
    multpoint: &mut bool,
    tol: f64,
) {
    let mut j = if *procf && *procl { slin.len() } else { 0 };
    while j < slin.len() {
        if slin[j].line_type == IntPatchIType::Analytic {
            let aligold_vertices = slin[j].a_curve.as_ref().map(|c| c.vertices.clone()).unwrap_or_default();
            let mut k = 0usize;
            while k < aligold_vertices.len() {
                let mut ptsol = aligold_vertices[k].clone();
                if !*procf {
                    let d = ptf.distance(ptsol.pnt.p);
                    if d <= tol {
                        ptsol.tolerance = tol;
                        if !ptsol.multiple {
                            *multpoint = true;
                            ptsol.multiple = true;
                        }
                        ptsol.param_on_line = first;
                        alig.vertices.push(ptsol.clone());
                        *procf = true;
                    }
                }
                if !*procl {
                    let d = ptl.distance(ptsol.pnt.p);
                    if d <= tol {
                        ptsol.tolerance = tol;
                        if !ptsol.multiple {
                            *multpoint = true;
                            ptsol.multiple = true;
                        }
                        ptsol.param_on_line = last;
                        alig.vertices.push(ptsol.clone());
                        *procl = true;
                    }
                }
                if *procf && *procl {
                    k = aligold_vertices.len();
                } else {
                    k += 1;
                }
            }
            if *procf && *procl {
                j = slin.len();
            } else {
                j += 1;
            }
        } else {
            j += 1;
        }
    }

    // Build a vertex for a 3D point on both surfaces.
    let make_vertex = |p: DVec3, param: f64, tol: f64| -> super::special_points::PatchPoint {
        let (u1, v1) = quad1.parameters(p);
        let (u2, v2) = quad2.parameters(p);
        super::special_points::PatchPoint {
            pnt: super::special_points::PntOn2S { p, u1, v1, u2, v2 },
            param_on_line: param,
            tolerance: tol,
            multiple: false,
            // OCCT IntPatch_Point::SetValue(Pt, Tol, Tangent) creates a point
            // "on no domain" (onS1 = onS2 = false); the ComputeVertexParameters
            // "remove first/last vertex not on any domain" step then drops the
            // endpoint vertex when it is not a surface-boundary crossing.
            on_dom_s1: false,
            on_dom_s2: false,
            arc_on_s1: None,
            arc_on_s2: None,
            param_on_arc1: 0.0,
            param_on_arc2: 0.0,
            is_vertex_on_s1: false,
            is_vertex_on_s2: false,
            transition_line_arc1: super::transitions::TypeTrans::Undecided,
            transition_line_arc2: super::transitions::TypeTrans::Undecided,
            transition_on_s1: super::transitions::TypeTrans::Undecided,
            transition_on_s2: super::transitions::TypeTrans::Undecided,
        }
    };

    let mut ptsol = make_vertex(ptf, 0.0, tol);
    if !*procf && !*procl {
        if ptf.distance(ptl) <= tol {
            ptsol.multiple = true;
            *multpoint = true;
            ptsol.param_on_line = first;
            alig.vertices.push(ptsol.clone());
            ptsol.param_on_line = last;
            alig.vertices.push(ptsol);
        } else {
            ptsol.param_on_line = first;
            alig.vertices.push(ptsol);
            ptsol = make_vertex(ptl, last, tol);
            alig.vertices.push(ptsol);
        }
    } else if !*procf {
        ptsol.param_on_line = first;
        alig.vertices.push(ptsol);
    } else if !*procl {
        ptsol = make_vertex(ptl, last, tol);
        alig.vertices.push(ptsol);
    }
}

const TOL_ANG: f64 = 1e-8;
