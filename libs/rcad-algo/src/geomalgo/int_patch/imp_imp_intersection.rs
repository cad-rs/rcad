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
            11 => self.int_pp(&q1, &q2, tol_tang),
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
                self.int_pcy(&q1, &q2, TOL_ANG, tol_tang, b_reverse, h);
                // OCCT L2592: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2596-2604: case 13/31 Plane/Cone
            13 | 31 => {
                self.int_pco(&q1, &q2, TOL_ANG, tol_tang, b_reverse);
                // OCCT L2602: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2606-2614: case 14/41 Plane/Sphere
            14 | 41 => {
                self.int_psp(&q1, &q2, TOL_ANG, tol_tang, b_reverse);
                // OCCT L2612: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2616-2624: case 15/51 Plane/Torus
            15 | 51 => {
                self.int_pto(&q1, &q2, tol_tang, b_reverse);
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
                self.int_cyco(&q1, &q2, s1, s2, tol_tang, b_reverse);
                // OCCT L2685: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2689-2697: case 24/42 Cylinder/Sphere
            24 | 42 => {
                self.int_cysp(&q1, &q2, s1, s2, tol_tang, b_reverse);
                // OCCT L2695: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2699-2707: case 25/52 Cylinder/Torus
            25 | 52 => {
                self.int_cyto(&q1, &q2, tol_tang, b_reverse);
                // OCCT L2705: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2709-2716: case 33 Cone/Cone
            33 => {
                self.int_coco(&q1, &q2, s1, s2, tol_tang);
                // OCCT L2714: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2718-2726: case 34/43 Cone/Sphere
            34 | 43 => {
                self.int_cosp(&q1, &q2, s1, s2, tol_tang, b_reverse);
                // OCCT L2724: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2728-2735: case 35/53 Cone/Torus
            35 | 53 => self.int_coto(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2737-2744: case 44 Sphere/Sphere
            44 => {
                self.int_spsp(&q1, &q2, tol_tang);
                // OCCT L2742: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2746-2754: case 45/54 Sphere/Torus
            45 | 54 => {
                self.int_spto(&q1, &q2, tol_tang, b_reverse);
                // OCCT L2752: bEmpty = empt
                b_empty = self.empt;
            }
            // OCCT L2756-2763: case 55 Torus/Torus
            55 => {
                self.int_toto(&q1, &q2, tol_tang);
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

        if !same_surf {
            // OCCT L2786-2787: AFunc.SetQuadric(quad2); AFunc.Set(S1);
            a_func.set_quadric(q2.clone());
            a_func.set_surface(s1.clone());
            // OCCT L2789: solrst.Perform(AFunc, D1, TolArc, TolTang);
            let mut d1 = super::so_on_bounds::Domain::new(uv1[0], uv1[1], uv1[2], uv1[3]);
            solrst.perform(&mut a_func, &mut d1, tol_arc, tol_tang, false);
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
            let mut d2 = super::so_on_bounds::Domain::new(uv2[0], uv2[1], uv2[2], uv2[3]);
            solrst.perform(&mut a_func, &mut d2, tol_arc, tol_tang, false);
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
                s1, s2, &pnt1, &mut self.slin, true, &q1, &q2, multpoint, tol_arc,
            );
            // OCCT L2912: PutPointsOnLine(S1, S2, pnt2, slin, false, D2, quad2, quad1, multpoint, TolArc);
            super::restriction::put_points_on_line(
                s1, s2, &pnt2, &mut self.slin, false, &q2, &q1, multpoint, tol_arc,
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

        // OCCT L2976-2995: ComputeVertexParameters(TolArc) for each GLine.
        // Sorts the vertices by parameter on line and removes coincident
        // duplicates (IntPatch_GLine::ComputeVertexParameters L421-).  RLine
        // vertices keep their insertion order (the RLine method is a no-op);
        // ALine lines are converted to WLine upstream.
        for line in self.slin.iter_mut() {
            let is_gline = matches!(
                line.line_type,
                IntPatchIType::Line
                    | IntPatchIType::Circle
                    | IntPatchIType::Ellipse
                    | IntPatchIType::Parabola
                    | IntPatchIType::Hyperbola
            );
            if !is_gline {
                continue;
            }
            line.vertices
                .sort_by(|a, b| a.param_on_line.partial_cmp(&b.param_on_line).unwrap_or(std::cmp::Ordering::Equal));
            let a_tol_pc = 1000.0 * rcad_kernel::precision::PCONFUSION;
            let mut dedup: Vec<IntPatchVertex> = Vec::with_capacity(line.vertices.len());
            for v in std::mem::take(&mut line.vertices) {
                if dedup.is_empty()
                    || (v.param_on_line - dedup.last().unwrap().param_on_line).abs() > a_tol_pc
                {
                    dedup.push(v);
                }
            }
            line.vertices = dedup;
        }

        // OCCT L2997-3040+: additional vertex placement for circles without vertices
    }

    // =====================================================================
    // Helper: run QuadQuadGeo and convert results to IntPatchLine
    // =====================================================================
    fn run_geo(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        // Determine which perform method to call based on surface types
        match (
            super::geom_abs_of_quadric(q1.surface_type()),
            super::geom_abs_of_quadric(q2.surface_type()),
        ) {
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Plane) => {
                geo.perform_plane_plane(q1, q2, 1e-8, tol)
            }
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cylinder)
            | (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Plane) => {
                geo.perform_plane_cylinder(q1, q2, 1e-8, tol, 0.0)
            }
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Sphere)
            | (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Plane) => {
                geo.perform_plane_sphere(q1, q2)
            }
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Plane) => {
                geo.perform_plane_cone(q1, q2, 1e-8, tol)
            }
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Plane) => {
                geo.perform_plane_torus(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Cylinder) => {
                geo.perform_cylinder_cylinder(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Sphere)
            | (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Cylinder) => {
                geo.perform_cylinder_sphere(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Cylinder) => {
                geo.perform_cylinder_cone(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Cylinder) => {
                geo.perform_cylinder_torus(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Sphere) => {
                geo.perform_sphere_sphere(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Sphere) => {
                geo.perform_sphere_cone(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Sphere) => {
                geo.perform_sphere_torus(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Cone) => {
                geo.perform_cone_cone(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Cone) => {
                geo.perform_cone_torus(q1, q2, tol)
            }
            (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Torus) => {
                geo.perform_torus_torus(q1, q2, tol)
            }
            _ => {}
        }
        if !geo.is_done() {
            return;
        }
        self.empt = geo.type_inter() == AnaResultType::Empty;
        self.tgte = false;
        self.oppo = false;
        // Convert curves
        for c in geo.to_curves() {
            let line_type = match &c {
                Curve3::Line(_) => IntPatchIType::Line,
                Curve3::Circle(_) => IntPatchIType::Circle,
                Curve3::Ellipse(_) => IntPatchIType::Ellipse,
                Curve3::Parabola(_) => IntPatchIType::Parabola,
                Curve3::Hyperbola(_) => IntPatchIType::Hyperbola,
                _ => IntPatchIType::Unknown,
            };
            let t_range = match &c {
                Curve3::Line(_) => [f64::NEG_INFINITY, f64::INFINITY],
                Curve3::Circle(_) | Curve3::Ellipse(_) => [0.0, std::f64::consts::TAU],
                _ => c.default_domain(),
            };
            self.slin.push(IntPatchLine {
                line_type,
                curve: c,
                t_range,
                pcurve1: None,
                pcurve2: None,
                tolerance: 1e-7,
                tang_tolerance: 1e-7,
                wline_pnts: Vec::new(),
                is_purging_allowed: false,
                wl_type: WLineType::Unknown,
                vertices: Vec::new(),
arc_on_s1: None,
                arc_on_s2: None,
                trans1: None,
                trans2: None,
                first_point: None,
                last_point: None,
                a_curve: None,
            });
        }
        self.my_done = IntStatus::OK;
    }

    // =====================================================================
    // 15 IntXX() functions — each calls QuadQuadGeo then post-processes
    // =====================================================================

    fn int_pp(&mut self, q1: &Quadric, q2: &Quadric, tol_tang: f64) {
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
            return;
        }
        self.my_done = IntStatus::OK;
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
    ) {
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
            return;
        }
        self.my_done = IntStatus::OK;
    }

    // OCCT L3446-3706: IntPCo — Plane/Cone intersection
    //       rcad IntPatchLine has no transition fields; transitions dropped.
    fn int_pco(&mut self, q1: &Quadric, q2: &Quadric, tol_ang: f64, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (plane, cone) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_plane_cone(plane, cone, tol_ang, tol);
        if !geo.is_done() {
            return;
        }
        let typint = geo.type_inter();
        let _nb_sol = geo.nb_solutions();
        self.empt = false;

        // OCCT L3489-3701: switch(typint)
        match typint {
            AnaResultType::Point => {
                // OCCT L3491-3501: IntAna_Point → spnt
                let psol = geo.point(1);
                let ptsol = IntPatchPoint {
                    p1: psol,
                    p2: psol,
                    u1: 0.0,
                    v1: 0.0,
                    u2: 0.0,
                    v2: 0.0,
                    tolerance: tol,
                };
                self.spnt.push(ptsol);
                // OCCT L3499: spnt.Append(ptsol)
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Line => {
                // OCCT L3503-3677: IntAna_Line
                //       Falls back to post_process_geo for line conversion.
                for c in geo.to_curves() {
                    let line_type = IntPatchIType::Line;
                    let t_range = [f64::NEG_INFINITY, f64::INFINITY];
                    self.slin.push(IntPatchLine {
                        line_type,
                        curve: c,
                        t_range,
                        vertices: Vec::new(),
arc_on_s1: None,
                        arc_on_s2: None,
                        trans1: None,
                        trans2: None,
                        first_point: None,
                        last_point: None,
                        a_curve: None,
                        pcurve1: None,
                        pcurve2: None,
                        tolerance: 1e-7,
                        tang_tolerance: 1e-7,
                        wline_pnts: Vec::new(),
                        is_purging_allowed: false,
                        wl_type: WLineType::Unknown,
                    });
                }
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Circle => {
                // OCCT L3680-3701: IntAna_Circle
                let c = geo.circle();
                self.slin.push(IntPatchLine {
                    line_type: IntPatchIType::Circle,
                    curve: Curve3::Circle(c),
                    t_range: [0.0, std::f64::consts::TAU],
                    vertices: Vec::new(),
arc_on_s1: None,
                    arc_on_s2: None,
                    trans1: None,
                    trans2: None,
                    first_point: None,
                    last_point: None,
                    a_curve: None,
                    pcurve1: None,
                    pcurve2: None,
                    tolerance: 1e-7,
                    tang_tolerance: 1e-7,
                    wline_pnts: Vec::new(),
                    is_purging_allowed: false,
                    wl_type: WLineType::Unknown,
                });
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Ellipse => {
                // OCCT L3704-3717: IntAna_Ellipse
                let e = geo.ellipse();
                self.slin.push(IntPatchLine {
                    line_type: IntPatchIType::Ellipse,
                    curve: Curve3::Ellipse(e),
                    t_range: [0.0, std::f64::consts::TAU],
                    vertices: Vec::new(),
arc_on_s1: None,
                    arc_on_s2: None,
                    trans1: None,
                    trans2: None,
                    first_point: None,
                    last_point: None,
                    a_curve: None,
                    pcurve1: None,
                    pcurve2: None,
                    tolerance: 1e-7,
                    tang_tolerance: 1e-7,
                    wline_pnts: Vec::new(),
                    is_purging_allowed: false,
                    wl_type: WLineType::Unknown,
                });
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Empty => {
                self.empt = true;
            }
            // OCCT L3340: default — Hyperbola, Parabola, etc.
            _ => {
                for c in geo.to_curves() {
                    let line_type = match &c {
                        Curve3::Line(_) => IntPatchIType::Line,
                        Curve3::Circle(_) => IntPatchIType::Circle,
                        Curve3::Ellipse(_) => IntPatchIType::Ellipse,
                        Curve3::Parabola(_) => IntPatchIType::Parabola,
                        Curve3::Hyperbola(_) => IntPatchIType::Hyperbola,
                        _ => IntPatchIType::Unknown,
                    };
                    let t_range = match &c {
                        Curve3::Circle(_) | Curve3::Ellipse(_) => [0.0, std::f64::consts::TAU],
                        Curve3::Line(_) => [f64::NEG_INFINITY, f64::INFINITY],
                        _ => [0.0, 1.0],
                    };
                    self.slin.push(IntPatchLine {
                        line_type,
                        curve: c,
                        t_range,
                        vertices: Vec::new(),
arc_on_s1: None,
                        arc_on_s2: None,
                        trans1: None,
                        trans2: None,
                        first_point: None,
                        last_point: None,
                        a_curve: None,
                        pcurve1: None,
                        pcurve2: None,
                        tolerance: 1e-7,
                        tang_tolerance: 1e-7,
                        wline_pnts: Vec::new(),
                        is_purging_allowed: false,
                        wl_type: WLineType::Unknown,
                    });
                }
                if !self.slin.is_empty() || !self.spnt.is_empty() {
                    self.empt = false;
                    self.my_done = IntStatus::OK;
                }
            }
        }
    }

    // OCCT L3352-3432: IntPSp — Plane/Sphere intersection
    fn int_psp(&mut self, q1: &Quadric, q2: &Quadric, _tol_ang: f64, _tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (plane, sphere) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_plane_sphere(plane, sphere);
        if !geo.is_done() {
            return;
        }
        let typint = geo.type_inter();
        self.empt = false;

        // OCCT L3391-3431: switch(typint)
        match typint {
            AnaResultType::Empty => {
                self.empt = true;
            }
            AnaResultType::Point => {
                // OCCT L3398-3407: IntAna_Point → spnt
                let psol = geo.point(1);
                let ptsol = IntPatchPoint {
                    p1: psol,
                    p2: psol,
                    u1: 0.0,
                    v1: 0.0,
                    u2: 0.0,
                    v2: 0.0,
                    tolerance: _tol,
                };
                self.spnt.push(ptsol);
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Circle => {
                // OCCT L3410-3431: IntAna_Circle → GLine
                let c = geo.circle();
                self.slin.push(IntPatchLine {
                    line_type: IntPatchIType::Circle,
                    curve: Curve3::Circle(c),
                    t_range: [0.0, std::f64::consts::TAU],
                    vertices: Vec::new(),
arc_on_s1: None,
                    arc_on_s2: None,
                    trans1: None,
                    trans2: None,
                    first_point: None,
                    last_point: None,
                    a_curve: None,
                    pcurve1: None,
                    pcurve2: None,
                    tolerance: 1e-7,
                    tang_tolerance: 1e-7,
                    wline_pnts: Vec::new(),
                    is_purging_allowed: false,
                    wl_type: WLineType::Unknown,
                });
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            _ => {
                for c in geo.to_curves() {
                    let line_type = match &c {
                        Curve3::Line(_) => IntPatchIType::Line,
                        Curve3::Circle(_) => IntPatchIType::Circle,
                        Curve3::Ellipse(_) => IntPatchIType::Ellipse,
                        Curve3::Parabola(_) => IntPatchIType::Parabola,
                        Curve3::Hyperbola(_) => IntPatchIType::Hyperbola,
                        _ => IntPatchIType::Unknown,
                    };
                    let t_range = match &c {
                        Curve3::Circle(_) | Curve3::Ellipse(_) => [0.0, std::f64::consts::TAU],
                        Curve3::Line(_) => [f64::NEG_INFINITY, f64::INFINITY],
                        _ => [0.0, 1.0],
                    };
                    self.slin.push(IntPatchLine {
                        line_type,
                        curve: c,
                        t_range,
                        vertices: Vec::new(),
arc_on_s1: None,
                        arc_on_s2: None,
                        trans1: None,
                        trans2: None,
                        first_point: None,
                        last_point: None,
                        a_curve: None,
                        pcurve1: None,
                        pcurve2: None,
                        tolerance: 1e-7,
                        tang_tolerance: 1e-7,
                        wline_pnts: Vec::new(),
                        is_purging_allowed: false,
                        wl_type: WLineType::Unknown,
                    });
                }
                if !self.slin.is_empty() || !self.spnt.is_empty() {
                    self.empt = false;
                    self.my_done = IntStatus::OK;
                }
            }
        }
    }

    fn int_pto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (plane, torus) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_plane_torus(plane, torus, tol);
        self.post_process_geo(&geo, tol);
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

    fn int_cyco(&mut self, q1: &Quadric, q2: &Quadric, s1: &Surface3, s2: &Surface3, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (cyl, cone) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_cone(cyl, cone, tol);
        if geo.is_done() && geo.type_inter() == AnaResultType::NoGeometricSolution
            && self.int_quad_quad_fallback(s1, s2, b_reverse)
        {
            return;
        }
        self.post_process_geo(&geo, tol);
    }

    fn int_cysp(&mut self, q1: &Quadric, q2: &Quadric, s1: &Surface3, s2: &Surface3, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        // perform_cylinder_sphere expects (cylinder, sphere)
        let (cyl, sph) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_sphere(cyl, sph, tol);
        if geo.is_done() && geo.type_inter() == AnaResultType::NoGeometricSolution
            && self.int_quad_quad_fallback(s1, s2, b_reverse)
        {
            return;
        }
        self.post_process_geo(&geo, tol);
    }

    fn int_cyto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (cyl, torus) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_torus(cyl, torus, tol);
        self.post_process_geo(&geo, tol);
    }

    fn int_coco(&mut self, q1: &Quadric, q2: &Quadric, s1: &Surface3, s2: &Surface3, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_cone_cone(q1, q2, tol);
        if geo.is_done() && geo.type_inter() == AnaResultType::NoGeometricSolution
            && self.int_quad_quad_fallback(s1, s2, false)
        {
            return;
        }
        self.post_process_geo(&geo, tol);
    }

    fn int_cosp(&mut self, q1: &Quadric, q2: &Quadric, s1: &Surface3, s2: &Surface3, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (sph, con) = if b_reverse { (q1, q2) } else { (q2, q1) };
        geo.perform_sphere_cone(sph, con, tol);
        if geo.is_done() && geo.type_inter() == AnaResultType::NoGeometricSolution
            && self.int_quad_quad_fallback(s1, s2, b_reverse)
        {
            return;
        }
        self.post_process_geo(&geo, tol);
    }

    /// OCCT IntCySp/IntCyCo/IntCoCo/IntCoSp L8263/L8465/L9022/L9349: when
    /// IntAna_QuadQuadGeo returns IntAna_NoGeometricSolution, fall back to
    /// IntAna_IntQuadQuad (the general cylinder/cone-quadric intersection).
    /// Returns true when the fallback ran (empty result is also a valid result).
    fn int_quad_quad_fallback(&mut self, s1: &Surface3, s2: &Surface3, b_reverse: bool) -> bool {
        // The explicit surface is the cylinder or cone.
        let (explicit, other) = if b_reverse { (s2, s1) } else { (s1, s2) };
        let Some(other_quad) = super::int_quad_quad::IntAnaQuadric::from_surface3(other) else {
            return false;
        };
        let mut iqq = super::int_quad_quad::IntQuadQuad::new();
        match explicit {
            Surface3::Cylinder(cyl) => iqq.perform_cylinder(cyl, &other_quad),
            Surface3::Cone(con) => iqq.perform_cone(con, &other_quad),
            _ => return false,
        }
        if !iqq.is_done() || iqq.nb_curves() == 0 {
            if iqq.is_done() && iqq.nb_curves() == 0 {
                // No curves is a valid result (Empty).
                self.empt = true;
                self.my_done = IntStatus::OK;
                return true;
            }
            return false;
        }
        // OCCT IntCySp/IntCyCo/IntCoCo/IntCoSp: wrap each IntAna_Curve result
        // in an IntPatch_ALine (ArcType = IntPatch_Analytic), add the endpoint
        // vertices via ProcessBounds, and append to slin.  The ALine is later
        // converted to a WLine by IntPatch_ALineToWLine::MakeWLine inside
        // IntPatch_Intersection::GeomGeomPerfom.
        let Some(q1) = Quadric::from_surface3(s1) else { return false; };
        let Some(q2) = Quadric::from_surface3(s2) else { return false; };
        for i in 0..iqq.nb_curves() {
            let Some(curve) = iqq.curve(i) else { continue };
            let mut curve = curve.clone();
            let d = curve.domain();
            let first = d[0];
            let last = d[1];
            let ptf = curve.value(first).unwrap_or(DVec3::ZERO);
            let ptl = curve.value(last).unwrap_or(DVec3::ZERO);
            let tol = 1e-7;
            // OCCT IntCySp ProcessBounds: add the endpoint vertices.
            let mut procf = false;
            let mut procl = false;
            let mut multpoint = false;
            process_bounds(
                &mut curve,
                &self.slin,
                &q1,
                &q2,
                &mut procf,
                ptf,
                first,
                &mut procl,
                ptl,
                last,
                &mut multpoint,
                tol,
            );
            self.slin.push(IntPatchLine {
                line_type: IntPatchIType::Analytic,
                curve: Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: ptf,
                    direction: if ptl != ptf {
                        (ptl - ptf).normalize_or_zero()
                    } else {
                        DVec3::X
                    },
                }),
                t_range: [first, last],
                pcurve1: None,
                pcurve2: None,
                tolerance: 1e-7,
                tang_tolerance: 1e-7,
                wline_pnts: Vec::new(),
                is_purging_allowed: false,
                wl_type: WLineType::ImpImp,
                vertices: Vec::new(),
arc_on_s1: None,
                arc_on_s2: None,
                trans1: None,
                trans2: None,
                first_point: None,
                last_point: None,
                a_curve: Some(curve),
            });
        }
        self.empt = false;
        self.my_done = IntStatus::OK;
        true
    }

    fn int_coto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, _b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_cone_torus(q1, q2, tol);
        self.post_process_geo(&geo, tol);
    }

    fn int_spsp(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_sphere_sphere(q1, q2, tol);
        if !geo.is_done() {
            return;
        }
        match geo.type_inter() {
            AnaResultType::Same => {
                                self.same_surf = true;
self.empt = false;
                self.tgte = true;
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Circle => {
                let c = geo.circle();
                self.slin.push(IntPatchLine {
                    line_type: IntPatchIType::Circle,
                    curve: Curve3::Circle(c),
                    t_range: [0.0, std::f64::consts::TAU],
                    pcurve1: None,
                    pcurve2: None,
                    tolerance: 1e-7,
                    tang_tolerance: 1e-7,
                    wline_pnts: Vec::new(),
                    is_purging_allowed: false,
                    wl_type: WLineType::Unknown,
                    vertices: Vec::new(),
arc_on_s1: None,
                    arc_on_s2: None,
                    trans1: None,
                    trans2: None,
                    first_point: None,
                    last_point: None,
                    a_curve: None,
                });
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            _ => {}
        }
    }

    fn int_spto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, _b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_sphere_torus(q1, q2, tol);
        self.post_process_geo(&geo, tol);
    }

    fn int_toto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_torus_torus(q1, q2, tol);
        self.post_process_geo(&geo, tol);
    }

    // Post-process QuadQuadGeo results into self.slin / self.spnt.
    // Used by IntXX functions where OCCT creates GLine with transitions
    // (rcad IntPatchLine drops transitions).
    fn post_process_geo(&mut self, geo: &QuadQuadGeo, tol: f64) {
        if !geo.is_done() {
            return;
        }
        let typ = geo.type_inter();
        // OCCT: Same → tgte=true
        if typ == AnaResultType::Same {
            self.empt = false;
            self.tgte = true;
            self.my_done = IntStatus::OK;
            return;
        }
        self.empt = typ == AnaResultType::Empty;
        self.tgte = false;
        self.oppo = false;
        if self.empt && geo.nb_solutions() == 0 {
            return;
        }
        // OCCT L678: Point → spnt
        if typ == AnaResultType::Point {
            let psol = geo.point(1);
            self.spnt.push(IntPatchPoint {
                p1: psol,
                p2: psol,
                u1: 0.0,
                v1: 0.0,
                u2: 0.0,
                v2: 0.0,
                tolerance: tol,
            });
            self.my_done = IntStatus::OK;
            return;
        }
        for c in geo.to_curves() {
            let line_type = match &c {
                Curve3::Line(_) => IntPatchIType::Line,
                Curve3::Circle(_) => IntPatchIType::Circle,
                Curve3::Ellipse(_) => IntPatchIType::Ellipse,
                Curve3::Parabola(_) => IntPatchIType::Parabola,
                Curve3::Hyperbola(_) => IntPatchIType::Hyperbola,
                _ => IntPatchIType::Unknown,
            };
            let t_range = match &c {
                Curve3::Circle(_) | Curve3::Ellipse(_) => [0.0, std::f64::consts::TAU],
                Curve3::Line(_) => [f64::NEG_INFINITY, f64::INFINITY],
                _ => [0.0, 1.0],
            };
            // OCCT L362-370: line created by QuadQuadGeo; vertices are added
            // later by PutPointsOnLine (intpatch_intersection post-process).
            self.slin.push(IntPatchLine {
                line_type,
                curve: c,
                t_range,
                vertices: Vec::new(),
arc_on_s1: None,
                arc_on_s2: None,
                trans1: None,
                trans2: None,
                first_point: None,
                last_point: None,
                a_curve: None,
                pcurve1: None,
                pcurve2: None,
                tolerance: 1e-7,
                tang_tolerance: 1e-7,
                wline_pnts: Vec::new(),
                is_purging_allowed: false,
                wl_type: WLineType::Unknown,
            });
        }
        self.my_done = IntStatus::OK;
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
            on_dom_s1: true,
            on_dom_s2: true,
            arc_on_s1: None,
            arc_on_s2: None,
            param_on_arc1: 0.0,
            param_on_arc2: 0.0,
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

/// Analytic UV inversion of a 3D point on a surface (ElSLib::Parameters).
fn surf_uv(surf: &Surface3, p: DVec3) -> (f64, f64) {
    use rcad_kernel::geom::SurfaceEval;
    match surf {
        Surface3::Plane(pl) => {
            let d = p - pl.origin;
            (d.dot(pl.u_dir), d.dot(pl.v_dir))
        }
        Surface3::Cylinder(c) => {
            let uv = c.world_to_uv(p);
            (uv.x, uv.y)
        }
        Surface3::Sphere(s) => {
            let uv = s.world_to_uv(p);
            (uv.x, uv.y)
        }
        Surface3::Cone(c) => {
            let uv = c.world_to_uv(p);
            (uv.x, uv.y)
        }
        Surface3::Torus(t) => {
            let uv = t.world_to_uv(p);
            (uv.x, uv.y)
        }
        _ => (0.0, 0.0),
    }
}
