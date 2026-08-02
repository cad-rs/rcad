// OCCT ElCLib — parameter/value/derivative for 3D conic curves (gp_Lin,
// gp_Circ, gp_Elips, gp_Parab, gp_Hypr), used by IntPatch_ImpImpIntersection
// IntXX functions and the restriction processing.
//
// 1:1 Rust translation of the ElCLib methods used there.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, CurveEval, Ellipse3, Hyperbola3, Line3, Parabola3};

/// OCCT ElCLib::Parameter(gp_Lin, gp_Pnt) — parameter of a point on a line.
pub fn line_parameter(line: &Line3, p: DVec3) -> f64 {
    (p - line.origin).dot(line.direction.normalize_or_zero())
}

/// OCCT ElCLib::Value(para, gp_Lin).
pub fn line_value(line: &Line3, t: f64) -> DVec3 {
    line.point_at(t)
}

/// OCCT ElCLib::D1(para, gp_Lin, P, V).
pub fn line_d1(line: &Line3, t: f64) -> (DVec3, DVec3) {
    (line.point_at(t), line.direction)
}

/// OCCT ElCLib::Parameter(gp_Circ, gp_Pnt) — angle parameter of a point on a circle.
pub fn circle_parameter(circ: &Circle3, p: DVec3) -> f64 {
    let d = p - circ.center;
    let x = d.dot(circ.x_dir.normalize_or_zero());
    let y = d.dot(circ.y_dir.normalize_or_zero());
    y.atan2(x)
}

/// OCCT ElCLib::Value(para, gp_Circ).
pub fn circle_value(circ: &Circle3, t: f64) -> DVec3 {
    circ.point_at(t)
}

/// OCCT ElCLib::D1(para, gp_Circ, P, V) — point and first derivative.
pub fn circle_d1(circ: &Circle3, t: f64) -> (DVec3, DVec3) {
    let x = circ.x_dir.normalize_or_zero();
    let y = circ.y_dir.normalize_or_zero();
    let p = circ.center + circ.radius * (t.cos() * x + t.sin() * y);
    let v = circ.radius * (-t.sin() * x + t.cos() * y);
    (p, v)
}

/// OCCT ElCLib::DN(para, gp_Circ, 1) — first derivative of the circle.
pub fn circle_dn(circ: &Circle3, t: f64) -> DVec3 {
    let x = circ.x_dir.normalize_or_zero();
    let y = circ.y_dir.normalize_or_zero();
    circ.radius * (-t.sin() * x + t.cos() * y)
}

/// OCCT ElCLib::Parameter(gp_Elips, gp_Pnt) — angle parameter of a point on an ellipse.
pub fn ellipse_parameter(e: &Ellipse3, p: DVec3) -> f64 {
    let d = p - e.center;
    let x = d.dot(e.major_dir.normalize_or_zero());
    let y = d.dot(e.normal.cross(e.major_dir).normalize_or_zero());
    y.atan2(x)
}

/// OCCT ElCLib::Value(para, gp_Elips).
pub fn ellipse_value(e: &Ellipse3, t: f64) -> DVec3 {
    e.point_at(t)
}

/// OCCT ElCLib::D1(para, gp_Elips, P, V).
pub fn ellipse_d1(e: &Ellipse3, t: f64) -> (DVec3, DVec3) {
    (e.point_at(t), e.derivative_at(t))
}

/// OCCT ElCLib::Parameter(gp_Parab, gp_Pnt) — approximate parameter.
pub fn parabola_parameter(p: &Parabola3, _pt: DVec3) -> f64 {
    p.default_domain()[0]
}

/// OCCT ElCLib::Value(para, gp_Parab).
pub fn parabola_value(p: &Parabola3, t: f64) -> DVec3 {
    p.point_at(t)
}

/// OCCT ElCLib::D1(para, gp_Parab, P, V).
pub fn parabola_d1(p: &Parabola3, t: f64) -> (DVec3, DVec3) {
    (p.point_at(t), p.derivative_at(t))
}

/// OCCT ElCLib::Parameter(gp_Hypr, gp_Pnt) — approximate parameter.
pub fn hyperbola_parameter(h: &Hyperbola3, _pt: DVec3) -> f64 {
    h.default_domain()[0]
}

/// OCCT ElCLib::Value(para, gp_Hypr).
pub fn hyperbola_value(h: &Hyperbola3, t: f64) -> DVec3 {
    h.point_at(t)
}

/// OCCT ElCLib::D1(para, gp_Hypr, P, V).
pub fn hyperbola_d1(h: &Hyperbola3, t: f64) -> (DVec3, DVec3) {
    (h.point_at(t), h.derivative_at(t))
}

/// OCCT ElCLib::LineParameter(gp_Ax1, gp_Pnt) — parameter of a point along an axis line.
pub fn line_parameter_of_axis(loc: DVec3, dir: DVec3, p: DVec3) -> f64 {
    (p - loc).dot(dir.normalize_or_zero())
}

/// OCCT ElCLib::Parameter(para, Curve3) — generic conic parameter dispatch.
pub fn conic_parameter(curve: &Curve3, p: DVec3) -> f64 {
    match curve {
        Curve3::Line(l) => line_parameter(l, p),
        Curve3::Circle(c) => circle_parameter(c, p),
        Curve3::Ellipse(e) => ellipse_parameter(e, p),
        Curve3::Parabola(pa) => parabola_parameter(pa, p),
        Curve3::Hyperbola(h) => hyperbola_parameter(h, p),
        _ => 0.0,
    }
}

/// OCCT ElCLib::Value(para, Curve3) — generic conic value dispatch.
pub fn conic_value(curve: &Curve3, t: f64) -> DVec3 {
    match curve {
        Curve3::Line(l) => line_value(l, t),
        Curve3::Circle(c) => circle_value(c, t),
        Curve3::Ellipse(e) => ellipse_value(e, t),
        Curve3::Parabola(pa) => parabola_value(pa, t),
        Curve3::Hyperbola(h) => hyperbola_value(h, t),
        _ => DVec3::ZERO,
    }
}

/// OCCT ElCLib::D1(para, Curve3) — generic conic D1 dispatch.
pub fn conic_d1(curve: &Curve3, t: f64) -> (DVec3, DVec3) {
    match curve {
        Curve3::Line(l) => line_d1(l, t),
        Curve3::Circle(c) => circle_d1(c, t),
        Curve3::Ellipse(e) => ellipse_d1(e, t),
        Curve3::Parabola(pa) => parabola_d1(pa, t),
        Curve3::Hyperbola(h) => hyperbola_d1(h, t),
        _ => (DVec3::ZERO, DVec3::ZERO),
    }
}

/// OCCT ElCLib::DN(para, Curve3, 1) — generic conic first-derivative dispatch.
pub fn conic_dn(curve: &Curve3, t: f64) -> DVec3 {
    match curve {
        Curve3::Line(l) => l.direction,
        Curve3::Circle(c) => circle_dn(c, t),
        Curve3::Ellipse(e) => e.derivative_at(t),
        Curve3::Parabola(pa) => pa.derivative_at(t),
        Curve3::Hyperbola(h) => h.derivative_at(t),
        _ => DVec3::ZERO,
    }
}

/// OCCT SeamPosition (IntPatch_ImpImpIntersection.cxx L3919-3925): rebuild the
/// circle's local frame with its location, the quadric's Z direction and X
/// direction.
pub fn seam_position(circ: &mut Circle3, a_pos_loc: DVec3, a_dz: DVec3, a_dx: DVec3) {
    let dz = a_dz.normalize_or_zero();
    let dx = a_dx.normalize_or_zero();
    let dy = dz.cross(dx).normalize_or_zero();
    circ.center = a_pos_loc;
    circ.normal = dz;
    circ.x_dir = dx;
    circ.y_dir = dy;
}

/// OCCT AdjustToSeam(const gp_Cone/Cylinder/Torus, gp_Circ) (L3863-3915):
/// reposition the circle on the quadric's seam frame.
pub fn adjust_to_seam_quadric(
    circ: &mut Circle3,
    quad_loc: DVec3,
    quad_dz: DVec3,
    quad_dx: DVec3,
) {
    let a_ploc = circ.center;
    seam_position(circ, a_ploc, quad_dz, quad_dx);
}

/// OCCT AdjustToSeam(const gp_Sphere, gp_Circ, aTolAng) (L3875-3891): only when
/// the circle axis is parallel to the sphere axis.
pub fn adjust_to_seam_sphere(
    circ: &mut Circle3,
    quad_loc: DVec3,
    quad_dz: DVec3,
    quad_dx: DVec3,
    tol_ang: f64,
) {
    let a_dir_c = circ.normal.normalize_or_zero();
    let a_dir_q = quad_dz.normalize_or_zero();
    if a_dir_c.dot(a_dir_q).abs() > 1.0 - tol_ang {
        let a_ploc = circ.center;
        seam_position(circ, a_ploc, a_dir_q, quad_dx);
    }
}
