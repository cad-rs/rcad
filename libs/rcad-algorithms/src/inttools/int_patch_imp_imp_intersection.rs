//! OCCT-aligned: IntPatch_ImpImpIntersection — intersection of two analytic surfaces.
//!
//! OCCT IntPatch_ImpImpIntersection.hxx / .cxx
//!
//! Handles all 15 pair combinations of Plane, Cylinder, Sphere, Cone, Torus
//! by converting to IntSurf_Quadric and dispatching to IntAna_QuadQuadGeo.

use rcad_kernel::geom::{Curve3, CurveEval, Surface3};
use super::int_surf_quadric::Quadric;
use super::int_ana_quad_quad_geo::{QuadQuadGeo, AnaResultType};
use super::int_patch_line::IntPatchLine;
use super::int_patch_line::WLineType;
use super::int_patch_point::IntPatchPoint;
use super::int_patch_type::IntPatchIType;
use super::geom_abs_surface_type::GeomAbsSurfaceType;

/// OCCT IntPatch_ImpImpIntersection.hxx L35-49: IntStatus
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntStatus { OK, InfiniteSectionCurve, Fail }

/// OCCT-aligned: IntPatch_ImpImpIntersection
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
            empt: true, tgte: false, oppo: false,
            spnt: Vec::new(), slin: Vec::new(),
        }
    }

    pub fn is_done(&self) -> bool { self.my_done == IntStatus::OK }
    pub fn status(&self) -> IntStatus { self.my_done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn tangent_faces(&self) -> bool { self.tgte }
    pub fn opposite_faces(&self) -> bool { self.oppo }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn line(&self, i: usize) -> &IntPatchLine { &self.slin[i] }
    pub fn slin_ref(&self) -> &[IntPatchLine] { &self.slin }
    pub fn nb_points(&self) -> usize { self.spnt.len() }
    pub fn point(&self, i: usize) -> &IntPatchPoint { &self.spnt[i] }

    // =====================================================================
    // OCCT L73-79: Perform(S1, D1, S2, D2, TolArc, TolTang, ...)
    // rcad: Surface3 instead of Adaptor3d_Surface. No TopolTool.
    // =====================================================================
    pub fn perform(&mut self, s1: &Surface3, s2: &Surface3, tol_arc: f64, tol_tang: f64) {
        self.my_done = IntStatus::Fail;
        self.spnt.clear();
        self.slin.clear();
        self.empt = true; self.tgte = false; self.oppo = false;

        // OCCT: SetQuad — convert Surface3 to Quadric + type index
        let Some(q1) = Quadric::from_surface3(s1) else { return; };
        let Some(q2) = Quadric::from_surface3(s2) else { return; };
        let i_t1 = quad_type_index(&q1);
        let i_t2 = quad_type_index(&q2);
        if i_t1 == 0 || i_t2 == 0 { return; }

        let b_reverse = i_t1 > i_t2;
        let i_tt = i_t1 * 10 + i_t2;

        // OCCT: dispatch based on iTT
        match i_tt {
            11 => self.int_pp(&q1, &q2, tol_tang),
            12 | 21 => self.int_pcy(&q1, &q2, TOL_ANG, tol_tang, b_reverse),
            13 | 31 => self.int_pco(&q1, &q2, TOL_ANG, tol_tang, b_reverse),
            14 | 41 => self.int_psp(&q1, &q2, TOL_ANG, tol_tang, b_reverse),
            15 | 51 => self.int_pto(&q1, &q2, tol_tang, b_reverse),
            22 => self.int_cycy(&q1, &q2, tol_tang),
            23 | 32 => self.int_cyco(&q1, &q2, tol_tang, b_reverse),
            24 | 42 => self.int_cysp(&q1, &q2, tol_tang, b_reverse),
            25 | 52 => self.int_cyto(&q1, &q2, tol_tang, b_reverse),
            33 => self.int_coco(&q1, &q2, tol_tang),
            34 | 43 => self.int_cosp(&q1, &q2, tol_tang, b_reverse),
            35 | 53 => self.int_coto(&q1, &q2, tol_tang, b_reverse),
            44 => self.int_spsp(&q1, &q2, tol_tang),
            45 | 54 => self.int_spto(&q1, &q2, tol_tang, b_reverse),
            55 => self.int_toto(&q1, &q2, tol_tang),
            _ => {}
        }

        // OCCT post-processing (skipped: PutPointsOnLine, ProcessSegments, ProcessRLine)
    }

    // =====================================================================
    // Helper: run QuadQuadGeo and convert results to IntPatchLine
    // =====================================================================
    fn run_geo(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        // Determine which perform method to call based on surface types
        match (q1.surface_type(), q2.surface_type()) {
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Plane) =>
                geo.perform_plane_plane(q1, q2, 1e-8, tol),
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cylinder)
            | (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Plane) =>
                geo.perform_plane_cylinder(q1, q2, 1e-8, tol, 0.0),
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Sphere)
            | (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Plane) =>
                geo.perform_plane_sphere(q1, q2),
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Plane) =>
                geo.perform_plane_cone(q1, q2, 1e-8, tol),
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Plane) =>
                geo.perform_plane_torus(q1, q2, tol),
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Cylinder) =>
                geo.perform_cylinder_cylinder(q1, q2, tol),
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Sphere)
            | (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Cylinder) =>
                geo.perform_cylinder_sphere(q1, q2, tol),
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Cylinder) =>
                geo.perform_cylinder_cone(q1, q2, tol),
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Cylinder) =>
                geo.perform_cylinder_torus(q1, q2, tol),
            (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Sphere) =>
                geo.perform_sphere_sphere(q1, q2, tol),
            (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Sphere) =>
                geo.perform_sphere_cone(q1, q2, tol),
            (GeomAbsSurfaceType::Sphere, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Sphere) =>
                geo.perform_sphere_torus(q1, q2, tol),
            (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Cone) =>
                geo.perform_cone_cone(q1, q2, tol),
            (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Torus)
            | (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Cone) =>
                geo.perform_cone_torus(q1, q2, tol),
            (GeomAbsSurfaceType::Torus, GeomAbsSurfaceType::Torus) =>
                geo.perform_torus_torus(q1, q2, tol),
            _ => {}
        }
        if !geo.is_done() { return; }
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
                line_type, curve: c, t_range,
                pcurve1: None, pcurve2: None,
                tolerance: 1e-7, tang_tolerance: 1e-7,
                wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: crate::inttools::int_patch_line::WLineType::Unknown,
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
        if !geo.is_done() { return; }
        match geo.type_inter() {
            AnaResultType::Same => { self.empt = false; self.tgte = true; self.my_done = IntStatus::OK; }
            AnaResultType::Line => {
                for c in geo.to_curves() {
                    self.slin.push(IntPatchLine {
                        line_type: IntPatchIType::Line, curve: c, t_range: [f64::NEG_INFINITY, f64::INFINITY],
                        pcurve1: None, pcurve2: None, tolerance: 1e-7, tang_tolerance: 1e-7,
                wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: crate::inttools::int_patch_line::WLineType::Unknown,
                    });
                }
                self.empt = false;
                self.my_done = IntStatus::OK;
            }
            _ => {}
        }
    }

    fn int_pcy(&mut self, q1: &Quadric, q2: &Quadric, _tol_ang: f64, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (plane, cyl) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_plane_cylinder(plane, cyl, 1e-8, tol, 0.0);
        self.post_process_geo(&geo);
    }

    fn int_pco(&mut self, q1: &Quadric, q2: &Quadric, tol_ang: f64, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (plane, cone) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_plane_cone(plane, cone, tol_ang, tol);
        self.post_process_geo(&geo);
    }

    fn int_psp(&mut self, q1: &Quadric, q2: &Quadric, _tol_ang: f64, _tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        // perform_plane_sphere expects (plane, sphere). b_reverse=true means q1=Sphere, q2=Plane.
        let (plane, sphere) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_plane_sphere(plane, sphere);
        self.post_process_geo(&geo);
    }

    fn int_pto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (plane, torus) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_plane_torus(plane, torus, tol);
        self.post_process_geo(&geo);
    }

    fn int_cycy(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_cylinder_cylinder(q1, q2, tol);
        if !geo.is_done() { return; }
        match geo.type_inter() {
            AnaResultType::Same => { self.empt = false; self.tgte = true; self.my_done = IntStatus::OK; }
            AnaResultType::Line => {
                for i in 1..=geo.nb_solutions() {
                    let l = geo.line(i);
                    self.slin.push(IntPatchLine {
                        line_type: IntPatchIType::Line,
                        curve: Curve3::Line(l), t_range: [-1e10, 1e10],
                        pcurve1: None, pcurve2: None,
                        tolerance: 1e-7, tang_tolerance: 1e-7,
                wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: crate::inttools::int_patch_line::WLineType::Unknown,
                    });
                }
                self.empt = false; self.my_done = IntStatus::OK;
            }
            AnaResultType::Ellipse => {
                let e = geo.ellipse();
                self.slin.push(IntPatchLine {
                    line_type: IntPatchIType::Ellipse,
                    curve: Curve3::Ellipse(e), t_range: [0.0, std::f64::consts::TAU],
                    pcurve1: None, pcurve2: None,
                    tolerance: 1e-7, tang_tolerance: 1e-7,
                wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: crate::inttools::int_patch_line::WLineType::Unknown,
                });
                self.empt = false; self.my_done = IntStatus::OK;
            }
            _ => {}
        }
    }

    fn int_cyco(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (cyl, cone) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_cone(cyl, cone, tol);
        self.post_process_geo(&geo);
    }

    fn int_cysp(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        // perform_cylinder_sphere expects (cylinder, sphere)
        let (cyl, sph) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_sphere(cyl, sph, tol);
        self.post_process_geo(&geo);
    }

    fn int_cyto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (cyl, torus) = if b_reverse { (q2, q1) } else { (q1, q2) };
        geo.perform_cylinder_torus(cyl, torus, tol);
        self.post_process_geo(&geo);
    }

    fn int_coco(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_cone_cone(q1, q2, tol);
        self.post_process_geo(&geo);
    }

    fn int_cosp(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        let (sph, con) = if b_reverse { (q1, q2) } else { (q2, q1) };
        geo.perform_sphere_cone(sph, con, tol);
        self.post_process_geo(&geo);
    }

    fn int_coto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, _b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_cone_torus(q1, q2, tol);
        self.post_process_geo(&geo);
    }

    fn int_spsp(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_sphere_sphere(q1, q2, tol);
        if !geo.is_done() { return; }
        match geo.type_inter() {
            AnaResultType::Same => { self.empt = false; self.tgte = true; self.my_done = IntStatus::OK; }
            AnaResultType::Circle => {
                let c = geo.circle();
                self.slin.push(IntPatchLine {
                    line_type: IntPatchIType::Circle,
                    curve: Curve3::Circle(c), t_range: [0.0, std::f64::consts::TAU],
                    pcurve1: None, pcurve2: None,
                    tolerance: 1e-7, tang_tolerance: 1e-7,
                wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: crate::inttools::int_patch_line::WLineType::Unknown,
                });
                self.empt = false; self.my_done = IntStatus::OK;
            }
            _ => {}
        }
    }

    fn int_spto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64, _b_reverse: bool) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_sphere_torus(q1, q2, tol);
        self.post_process_geo(&geo);
    }

    fn int_toto(&mut self, q1: &Quadric, q2: &Quadric, tol: f64) {
        let mut geo = QuadQuadGeo::new();
        geo.perform_torus_torus(q1, q2, tol);
        self.post_process_geo(&geo);
    }

    // =====================================================================
    // Post-process QuadQuadGeo results into self.slin / self.spnt
    // =====================================================================
    fn post_process_geo(&mut self, geo: &QuadQuadGeo) {
        if !geo.is_done() { return; }
        self.empt = geo.type_inter() == AnaResultType::Empty;
        self.tgte = false;
        self.oppo = false;
        if self.empt && geo.nb_solutions() == 0 { return; }
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
                line_type, curve: c, t_range,
                pcurve1: None, pcurve2: None,
                tolerance: 1e-7, tang_tolerance: 1e-7,
                wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: crate::inttools::int_patch_line::WLineType::Unknown,
            });
        }
        self.my_done = IntStatus::OK;
    }
}

/// OCCT: SetQuad equivalent — return type index (1=Plane, 2=Cyl, 3=Cone, 4=Sphere, 5=Torus)
fn quad_type_index(q: &Quadric) -> i32 {
    match q.surface_type() {
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
