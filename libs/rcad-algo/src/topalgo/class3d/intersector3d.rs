// OCCT BRepClass3d_Intersector3d (BRepClass3d_Intersector3d.hxx / .cxx)
// 3D intersection of a line segment with a face.

use glam::DVec3;
use crate::topalgo::int_curve_surface::intersector::TransitionOnCurve;

/// OCCT BRepClass3d_Intersector3d — intersects a line segment with a face.
pub struct Intersector3d {
    u: f64,
    v: f64,
    w: f64,
    transition: TransitionOnCurve,
    done: bool,
    has_a_point: bool,
    state: u8,  // 0=UNKNOWN, 1=IN, 2=ON, 3=OUT
    point: DVec3,
    face_idx: usize,
}

impl Intersector3d {
    pub fn new() -> Self {
        Intersector3d {
            u: 0.0, v: 0.0, w: 0.0,
            transition: TransitionOnCurve::In,
            done: false, has_a_point: false,
            state: 0, point: DVec3::ZERO, face_idx: usize::MAX,
        }
    }

    /// OCCT: Perform(L, Prm, Tol, F) — intersect line L with face F.
    ///
    /// OCCT L29-146: uses IntCurveSurface_HInter + BRepClass_FaceClassifier.
    /// rcad: delegates to IntCurvesFace_Intersector (topalgo::int_curve_surface)
    /// which uses rcad-kernel::IntCS for curve-surface intersection.
    pub fn perform(&mut self, line_origin: DVec3, line_dir: DVec3,
                   _prm: f64, _tol: f64,
                   face_surface: &rcad_kernel::geom::Surface3,
                   face_idx: usize) {
        // Build a line curve for intersection
        let line_curve = rcad_kernel::geom::Curve3::Line(
            rcad_kernel::geom::Line3 { origin: line_origin, direction: line_dir });

        // OCCT: IntCurvesFace_Intersector intersects curve with face surface
        // rcad: use IntCurvesFace_Intersector via topalgo module
        let mut cf_intersector = crate::topalgo::int_curve_surface::intersector::Intersector::new();
        cf_intersector.perform(
            &line_curve, face_surface,
            0.0, 1.0, 0.0, 1.0, // default UV bounds
            false, false, 0.0, 0.0);

        self.done = cf_intersector.is_done();
        self.has_a_point = cf_intersector.has_a_point();
        self.u = cf_intersector.u_parameter();
        self.v = cf_intersector.v_parameter();
        self.w = cf_intersector.w_parameter();
        self.point = cf_intersector.pnt();
        self.transition = cf_intersector.transition();
        self.state = cf_intersector.state();
        self.face_idx = face_idx;
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn has_a_point(&self) -> bool { self.has_a_point }
    pub fn u_parameter(&self) -> f64 { self.u }
    pub fn v_parameter(&self) -> f64 { self.v }
    pub fn w_parameter(&self) -> f64 { self.w }
    pub fn pnt(&self) -> DVec3 { self.point }
    pub fn transition(&self) -> TransitionOnCurve { self.transition }
    pub fn state(&self) -> u8 { self.state }
    pub fn face_idx(&self) -> usize { self.face_idx }
}
