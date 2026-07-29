// OCCT BRepIntCurveSurface_Inter (BRepIntCurveSurface_Inter.hxx / .cxx)
// Intersection between a face/Shape and a curve.
// Provides iteration over intersection points.

use rcad_kernel::geom::{Curve3, CurveEval};

/// OCCT IntCurveSurface_TransitionOnCurve — transition type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransitionOnCurve {
    In, Out, Tangent,
}

/// OCCT IntCurveSurface_IntersectionPoint (IntCurveSurface_IntersectionPoint.hxx).
#[derive(Debug, Clone)]
pub struct IntersectionPoint {
    pub point: glam::DVec3,
    pub u: f64,
    pub v: f64,
    pub w: f64,
    pub transition: TransitionOnCurve,
}

/// OCCT BRepIntCurveSurface_Inter — intersects a curve with a face/shell/solid.
///
/// rcad: simplified — intersects a curve with a single face surface.
pub struct Inter {
    // OCCT fields
    my_cur_curve: Option<Curve3>,
    my_tolerance: f64,
    // Results
    current_index: usize,
    current_points: Vec<IntersectionPoint>,
    // State
    current_state: u8,  // 0=UNKNOWN, 1=IN, 2=ON, 3=OUT
    current_u: f64,
    current_v: f64,
    current_w: f64,
    current_point: glam::DVec3,
    current_transition: TransitionOnCurve,
}

impl Inter {
    pub fn new() -> Self {
        Inter {
            my_cur_curve: None, my_tolerance: 1e-7,
            current_index: 0, current_points: Vec::new(),
            current_state: 0, current_u: 0.0, current_v: 0.0,
            current_w: 0.0, current_point: glam::DVec3::ZERO,
            current_transition: TransitionOnCurve::In,
        }
    }

    /// OCCT: Load(Shape, Tol) — load the shape and tolerance.
    pub fn load(&mut self, _the_shape: &rcad_kernel::topo_shape::Shape, the_tol: f64) {
        self.clear();
        self.my_tolerance = the_tol;
        // rcad: store shape for face iteration (simplified — single face)
    }

    /// OCCT: Init(Curve) — initialize with curve, find all intersections.
    pub fn init_curve(&mut self, curve: &Curve3, face_surface: &rcad_kernel::geom::Surface3,
                      u_min: f64, u_max: f64, v_min: f64, v_max: f64) {
        self.clear();
        self.my_cur_curve = Some(curve.clone());

        // OCCT: IntCurveSurface_HInter for curve-surface intersection
        let mut hics = rcad_kernel::base::geom_api::int_cs::IntCS::new();
        hics.perform(curve, face_surface);

        if hics.is_done() {
            let nb_pts = hics.nb_points();
            for idx in (1..=nb_pts).rev() {
                let pt = hics.point(idx);
                let in_face = pt.u >= u_min && pt.u <= u_max
                           && pt.v >= v_min && pt.v <= v_max;
                if in_face {
                    // OCCT L111-139: classify 2D point against face
                    // rcad: accept if within UV bounds
                    self.current_points.push(IntersectionPoint {
                        point: pt.point,
                        u: pt.u, v: pt.v, w: pt.w,
                        transition: TransitionOnCurve::In,
                    });
                }
            }
        }
    }

    /// OCCT: More() — returns true if there are more intersection points.
    pub fn more(&self) -> bool {
        self.current_index < self.current_points.len()
    }

    /// OCCT: Next() — advance to next intersection point.
    pub fn next(&mut self) {
        if self.more() {
            let pt = &self.current_points[self.current_index];
            self.current_u = pt.u;
            self.current_v = pt.v;
            self.current_w = pt.w;
            self.current_point = pt.point;
            self.current_transition = pt.transition;
            self.current_state = 1; // IN
            self.current_index += 1;
        }
    }

    /// OCCT: Point() — current intersection point data.
    pub fn point(&self) -> &IntersectionPoint {
        &self.current_points[self.current_index.min(1).max(1) - 1]
    }

    pub fn clear(&mut self) {
        self.current_index = 0;
        self.current_points.clear();
        self.current_state = 0;
        self.current_u = 0.0;
        self.current_v = 0.0;
    }

    // Accessors matching OCCT interface
    pub fn current_u(&self) -> f64 { self.current_u }
    pub fn current_v(&self) -> f64 { self.current_v }
    pub fn current_w(&self) -> f64 { self.current_w }
    pub fn current_point(&self) -> glam::DVec3 { self.current_point }
    pub fn current_transition(&self) -> TransitionOnCurve { self.current_transition }
    pub fn current_state(&self) -> u8 { self.current_state }
}
