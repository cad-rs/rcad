//! IntPatch_ImpImpIntersection — intersection of two analytic surfaces.
//!
//! OCCT IntPatch_ImpImpIntersection.hxx / .cxx
//!
//! Handles all 15 pair combinations of Plane, Cylinder, Sphere, Cone, Torus
//! by converting to IntSurf_Quadric and dispatching to IntAna_QuadQuadGeo.

use super::GeomAbsSurfaceType;
use super::{AnaResultType, IntPatchLine, IntPatchPoint, IntPatchIType, QuadQuadGeo, WLineType};
use crate::topalgo::int_surf::quadric::Quadric;
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
    // rcad: Surface3 instead of Adaptor3d_Surface. No TopolTool.
    // =====================================================================
    pub fn perform(&mut self, s1: &Surface3, s2: &Surface3, tol_arc: f64, tol_tang: f64) {
        // OCCT L2525-2533: myDone=Fail, spnt/slin clear, empt/tgte/oppo init
        self.my_done = IntStatus::Fail;
        self.spnt.clear();
        self.slin.clear();
        self.empt = true;
        self.tgte = false;
        self.oppo = false;

        // OCCT L2529: isPostProcessingRequired = true
        // OCCT L2535-2546: all1, all2, SameSurf, multpoint, nosolonS1, nosolonS2
        //   edg1, edg2, pnt1, pnt2 — IntPatch_TheSOnBounds (needs TopolTool)
        let _is_post_processing_required = true;

        // OCCT L2548-2556: SetQuad — convert Surface3 to Quadric + type index
        let Some(q1) = Quadric::from_surface3(s1) else {
            return;
        };
        let Some(q2) = Quadric::from_surface3(s2) else {
            return;
        };
        // OCCT L2555: typs1, typs2 — surface type enums (rcad: inferred from Quadric)
        let _typs1 = q1.surface_type();
        let _typs2 = q2.surface_type();
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
            12 | 21 => self.int_pcy(&q1, &q2, TOL_ANG, tol_tang, b_reverse),
            // OCCT L2596-2604: case 13/31 Plane/Cone
            13 | 31 => self.int_pco(&q1, &q2, TOL_ANG, tol_tang, b_reverse),
            // OCCT L2606-2614: case 14/41 Plane/Sphere
            14 | 41 => self.int_psp(&q1, &q2, TOL_ANG, tol_tang, b_reverse),
            // OCCT L2616-2624: case 15/51 Plane/Torus
            15 | 51 => self.int_pto(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2626-2677: case 22 Cylinder/Cylinder (aBox1,aBox2,a2DTol)
            22 => self.int_cycy(&q1, &q2, tol_tang),
            // OCCT L2679-2687: case 23/32 Cylinder/Cone
            23 | 32 => self.int_cyco(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2689-2697: case 24/42 Cylinder/Sphere
            24 | 42 => self.int_cysp(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2699-2707: case 25/52 Cylinder/Torus
            25 | 52 => self.int_cyto(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2709-2716: case 33 Cone/Cone
            33 => self.int_coco(&q1, &q2, tol_tang),
            // OCCT L2718-2726: case 34/43 Cone/Sphere
            34 | 43 => self.int_cosp(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2728-2735: case 35/53 Cone/Torus
            35 | 53 => self.int_coto(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2737-2744: case 44 Sphere/Sphere
            44 => self.int_spsp(&q1, &q2, tol_tang),
            // OCCT L2746-2754: case 45/54 Sphere/Torus
            45 | 54 => self.int_spto(&q1, &q2, tol_tang, b_reverse),
            // OCCT L2756-2763: case 55 Torus/Torus
            55 => self.int_toto(&q1, &q2, tol_tang),
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
        // OCCT: solrst.Perform(AFunc, D1/D2, TolArc, TolTang) — boundary intersection
        //       PutPointsOnLine, ProcessSegments, ProcessRLine
        //       Boundary clipping is done upstream in IntPatch_Intersection.

        // OCCT L2936-2995: ComputeVertexParameters for each line
        for i in 0..self.slin.len() {
            // OCCT L2976-2995: ComputeVertexParameters(TolArc) for GLine/ALine/RLine
            //       Vertex parameters are computed in MakeCurve (upstream).
            let _ = i;
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
            });
        }
        self.my_done = IntStatus::OK;
    }

    // =====================================================================
    // 15 IntXX() functions — each calls QuadQuadGeo then post-processes
    // =====================================================================

    fn int_pp(&mut self, q1: &Quadric, q2: &Quadric, tol_tang: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_plane_plane(q1, q2, 1e-8, tol_tang);
        if !geo.is_done() {
            return;
        }
        match geo.type_inter() {
            AnaResultType::Same => {
                self.empt = false;
                self.tgte = true;
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Line => {
                for c in geo.to_curves() {
                    self.slin.push(IntPatchLine {
                        line_type: IntPatchIType::Line,
                        curve: c,
                        t_range: [f64::NEG_INFINITY, f64::INFINITY],
                        pcurve1: None,
                        pcurve2: None,
                        tolerance: 1e-7,
                        tang_tolerance: 1e-7,
                        wline_pnts: Vec::new(),
                        is_purging_allowed: false,
                        wl_type: WLineType::Unknown,
                        vertices: Vec::new(),
                    });
                }
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            _ => {}
        }
    }

    // OCCT L3157-3345: IntPCy — Plane/Cylinder intersection
    fn int_pcy(&mut self, q1: &Quadric, q2: &Quadric, tol_ang: f64, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (plane, cyl) = if b_reverse { (q2, q1) } else { (q1, q2) };
        // OCCT L3184: inter.Perform(Pl, Cy, Tolang, TolTang, H)
        geo.perform_plane_cylinder(plane, cyl, tol_ang, tol, 0.0);
        // OCCT L3185-3188: if (!inter.IsDone()) return false
        if !geo.is_done() {
            return;
        }
        // OCCT L3189-3191: typint, NbSol, Empty=false
        let typint = geo.type_inter();
        let _nb_sol = geo.nb_solutions();
        self.empt = false;

        // OCCT L3193-3343: switch(typint)
        match typint {
            // OCCT L3195-3198: case IntAna_Empty
            AnaResultType::Empty => {
                self.empt = true;
            }
            // OCCT L3200-3290: case IntAna_Line — 1 or 2 lines
            AnaResultType::Line => {
                // OCCT L3201: linsol = inter.Line(1)
                let linsol = geo.line(1);
                // OCCT L3258-3289: 2 lines (NbSol==2 path, no transition drop)
                self.slin.push(IntPatchLine::analytic(
                    IntPatchIType::Line,
                    Curve3::Line(linsol),
                    [f64::NEG_INFINITY, f64::INFINITY],
                ));
                if _nb_sol >= 2 {
                    let linsol2 = geo.line(2);
                    self.slin.push(IntPatchLine::analytic(
                        IntPatchIType::Line,
                        Curve3::Line(linsol2),
                        [f64::NEG_INFINITY, f64::INFINITY],
                    ));
                }
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            // OCCT L3293-3316: case IntAna_Circle
            AnaResultType::Circle => {
                // OCCT L3298: cirsol = inter.Circle(1)
                let cirsol = geo.circle();
                // GLine(cirsol, false, trans1, trans2) — transitions not in IntPatchLine
                self.slin.push(IntPatchLine::analytic(
                    IntPatchIType::Circle,
                    Curve3::Circle(cirsol),
                    [0.0, std::f64::consts::TAU],
                ));
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            // OCCT L3319-3337: case IntAna_Ellipse
            AnaResultType::Ellipse => {
                let elipsol = geo.ellipse();
                self.slin.push(IntPatchLine::analytic(
                    IntPatchIType::Ellipse,
                    Curve3::Ellipse(elipsol),
                    [0.0, std::f64::consts::TAU],
                ));
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            // OCCT L3340-3342: default
            _ => {}
        }
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

    fn int_cycy(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_cylinder_cylinder(q1, q2, tol);
        if !geo.is_done() {
            return;
        }
        match geo.type_inter() {
            AnaResultType::Same => {
                self.empt = false;
                self.tgte = true;
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Line => {
                for i in 1..=geo.nb_solutions() {
                    let l = geo.line(i);
                    self.slin.push(IntPatchLine {
                        line_type: IntPatchIType::Line,
                        curve: Curve3::Line(l),
                        t_range: [-1e10, 1e10],
                        pcurve1: None,
                        pcurve2: None,
                        tolerance: 1e-7,
                        tang_tolerance: 1e-7,
                        wline_pnts: Vec::new(),
                        is_purging_allowed: false,
                        wl_type: WLineType::Unknown,
                        vertices: Vec::new(),
                    });
                }
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            AnaResultType::Ellipse => {
                let e = geo.ellipse();
                self.slin.push(IntPatchLine {
                    line_type: IntPatchIType::Ellipse,
                    curve: Curve3::Ellipse(e),
                    t_range: [0.0, std::f64::consts::TAU],
                    pcurve1: None,
                    pcurve2: None,
                    tolerance: 1e-7,
                    tang_tolerance: 1e-7,
                    wline_pnts: Vec::new(),
                    is_purging_allowed: false,
                    wl_type: WLineType::Unknown,
                    vertices: Vec::new(),
                });
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            _ => {}
        }
    }

    fn int_cyco(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (cyl, cone) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_cone(cyl, cone, tol);
        self.post_process_geo(&geo, tol);
    }

    fn int_cysp(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        // perform_cylinder_sphere expects (cylinder, sphere)
        let (cyl, sph) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_sphere(cyl, sph, tol);
        self.post_process_geo(&geo, tol);
    }

    fn int_cyto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (cyl, torus) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_torus(cyl, torus, tol);
        self.post_process_geo(&geo, tol);
    }

    fn int_coco(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_cone_cone(q1, q2, tol);
        self.post_process_geo(&geo, tol);
    }

    fn int_cosp(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (sph, con) = if b_reverse { (q1, q2) } else { (q2, q1) };
        geo.perform_sphere_cone(sph, con, tol);
        self.post_process_geo(&geo, tol);
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
const TOL_ANG: f64 = 1e-8;
