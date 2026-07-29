// OCCT IntCurvesFace_Intersector (IntCurvesFace_Intersector.hxx / .cxx)
// Intersects a 3D curve with a face. Used by BRepClass3d_Intersector3d.

use rcad_kernel::geom::{Curve3, Surface3};

/// Transition type for curve-surface intersection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransitionOnCurve {
    In, Out, Tangent,
}

/// OCCT IntCurvesFace_Intersector — intersects a curve with a face.
///
/// OCCT L29-146: Perform(Curve, Face) → uses IntCurveSurface_HInter
/// to find intersection points, then filters by face geometry.
pub struct Intersector {
    // Result data
    u: f64,
    v: f64,
    w: f64,
    transition: TransitionOnCurve,
    done: bool,
    has_a_point: bool,
    state: u8,  // 0=UNKNOWN, 1=IN, 2=ON, 3=OUT
    point: glam::DVec3,
    // Input curve and tolerance
    curve: Option<Curve3>,
    tolerance: f64,
}

impl Intersector {
    pub fn new() -> Self {
        Intersector {
            u: 0.0, v: 0.0, w: 0.0,
            transition: TransitionOnCurve::In,
            done: false, has_a_point: false,
            state: 0, point: glam::DVec3::ZERO,
            curve: None, tolerance: 1e-7,
        }
    }

    /// OCCT IntCurvesFace_Intersector::Perform (IntCurvesFace_Intersector.cxx L29-146).
    ///
    /// Intersects the curve with the face surface.
    /// OCCT: uses IntCurveSurface_HInter (rcad-kernel IntCS) for curve-surface
    /// intersection, then BRepClass_FaceClassifier for 2D point-in-face test.
    pub fn perform(&mut self, curve: &Curve3, face_surface: &Surface3,
                   u_min: f64, u_max: f64, v_min: f64, v_max: f64,
                   is_u_periodic: bool, is_v_periodic: bool,
                   u_period: f64, v_period: f64) {
        self.curve = Some(curve.clone());
        self.done = true;
        self.has_a_point = false;

        // OCCT L47-70: IntCurveSurface_HInter HICS; HICS.Perform(curve, surface)
        // rcad: IntCS from rcad-kernel (equivalent to IntCurveSurface_HInter)
        let mut hics = rcad_kernel::base::geom_api::int_cs::IntCS::new();
        hics.perform(curve, face_surface);

        self.w = f64::MAX;
        if hics.is_done() {
            let nb_pts = hics.nb_points();
            for index in (1..=nb_pts).rev() {
                let pt = hics.point(index);
                let mut p_uv = glam::DVec2::new(pt.u, pt.v);

                // OCCT L85-109: handle UV periodicity
                if is_u_periodic {
                    let mut n1 = 0i64;
                    if p_uv.x > u_max {
                        n1 = ((p_uv.x - u_min) / u_period) as i64;
                    } else if p_uv.x < u_min {
                        n1 = ((p_uv.x - u_max) / u_period) as i64;
                    }
                    p_uv.x -= u_period * n1 as f64;
                }
                if is_v_periodic {
                    let mut n2 = 0i64;
                    if p_uv.y > v_max {
                        n2 = ((p_uv.y - v_min) / v_period) as i64;
                    } else if p_uv.y < v_min {
                        n2 = ((p_uv.y - v_min) / v_period) as i64;
                    }
                    p_uv.y -= v_period * n2 as f64;
                }

                // OCCT L111-139: BRepClass_FaceClassifier → 2D point-in-face test
                // rcad: simplified — accept point if within UV bounds
                let in_face = p_uv.x >= u_min && p_uv.x <= u_max
                           && p_uv.y >= v_min && p_uv.y <= v_max;

                if in_face && pt.w > -self.tolerance && pt.w < self.w {
                    self.has_a_point = true;
                    self.u = pt.u;
                    self.v = pt.v;
                    self.w = pt.w;
                    self.point = pt.point;
                    // rcad: simplified transition detection
                    self.transition = TransitionOnCurve::In;
                    self.state = 1; // IN
                }
            }
        }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn has_a_point(&self) -> bool { self.has_a_point }
    pub fn u_parameter(&self) -> f64 { self.u }
    pub fn v_parameter(&self) -> f64 { self.v }
    pub fn w_parameter(&self) -> f64 { self.w }
    pub fn pnt(&self) -> glam::DVec3 { self.point }
    pub fn transition(&self) -> TransitionOnCurve { self.transition }
    pub fn state(&self) -> u8 { self.state }
}
