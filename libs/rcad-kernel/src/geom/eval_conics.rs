//! Per-kind evaluation impls for the analytic curve kinds: the OCCT
//! `Geom_Conic` / `Geom2d_Conic` group (line, circle, ellipse, hyperbola,
//! parabola), the helix and sine wave, and the 2D analytic PCurve kinds
//! (involute, spirals, sine wave).
//!
//! Extracted verbatim from `geom/eval.rs` (project Rule 5 file-size split);
//! the evaluation traits stay re-exported from `geom/eval.rs`.
use crate::geom::*;
use crate::core::precision::INFINITE_VALUE;
use std::f64::consts::PI;

impl ConicEval for Circle3 {
    fn position(&self) -> DVec3 { self.center }
    fn normal(&self) -> DVec3 { self.normal }
    fn eccentricity(&self) -> f64 { 0.0 }
    fn x_axis(&self) -> DVec3 { self.x_dir }
    fn y_axis(&self) -> DVec3 { self.y_dir }
}

impl ConicEval for Ellipse3 {
    fn position(&self) -> DVec3 { self.center }
    fn normal(&self) -> DVec3 { self.normal }
    fn eccentricity(&self) -> f64 {
        let a = self.major_radius;
        let b = self.minor_radius;
        if a.abs() < 1e-15 { 1.0 } else { (1.0 - (b * b) / (a * a)).sqrt() }
    }
    fn x_axis(&self) -> DVec3 { self.major_dir }
    fn y_axis(&self) -> DVec3 { self.normal.cross(self.major_dir).normalize_or_zero() }
}

impl ConicEval for Hyperbola3 {
    fn position(&self) -> DVec3 { self.center }
    fn normal(&self) -> DVec3 { self.normal }
    fn eccentricity(&self) -> f64 {
        let a = self.semi_major;
        let b = self.semi_minor;
        if a.abs() < 1e-15 { 1.0 } else { (1.0 + (b * b) / (a * a)).sqrt() }
    }
    fn x_axis(&self) -> DVec3 { self.major_dir }
    fn y_axis(&self) -> DVec3 { self.normal.cross(self.major_dir).normalize_or_zero() }
}

impl ConicEval for Parabola3 {
    fn position(&self) -> DVec3 { self.vertex }
    fn normal(&self) -> DVec3 { self.normal }
    fn eccentricity(&self) -> f64 { 1.0 }
    fn x_axis(&self) -> DVec3 { self.axis_dir }
    fn y_axis(&self) -> DVec3 { self.axis_dir.cross(self.normal).normalize_or_zero() }
}

/// OCCT-aligned: `Geom2d_Conic` intermediate abstract class.
///
/// Groups 2D conic curves (Circle, Ellipse, Hyperbola, Parabola).
// --- Conic2dEval implementations ---

impl Conic2dEval for Circle2d {
    fn position(&self) -> DVec2 { self.center }
    fn eccentricity(&self) -> f64 { 0.0 }
    fn x_axis(&self) -> DVec2 { self.x_dir }
    fn y_axis(&self) -> DVec2 { self.y_dir }
}

impl Conic2dEval for Ellipse2d {
    fn position(&self) -> DVec2 { self.center }
    fn eccentricity(&self) -> f64 {
        let a = self.major_radius;
        let b = self.minor_radius;
        if a.abs() < 1e-15 { 1.0 } else { (1.0 - (b * b) / (a * a)).sqrt() }
    }
    fn x_axis(&self) -> DVec2 { self.major_dir }
    fn y_axis(&self) -> DVec2 { self.minor_dir }
}

impl Conic2dEval for Parabola2d {
    fn position(&self) -> DVec2 { self.origin }
    fn eccentricity(&self) -> f64 { 1.0 }
    fn x_axis(&self) -> DVec2 { self.axis_dir }
    fn y_axis(&self) -> DVec2 { DVec2::new(-self.axis_dir.y, self.axis_dir.x) }
}

impl Conic2dEval for Hyperbola2d {
    fn position(&self) -> DVec2 { self.center }
    fn eccentricity(&self) -> f64 {
        let a = self.semi_major;
        let b = self.semi_minor;
        if a.abs() < 1e-15 { 1.0 } else { (1.0 + (b * b) / (a * a)).sqrt() }
    }
    fn x_axis(&self) -> DVec2 { self.major_dir }
    fn y_axis(&self) -> DVec2 { turn_2d(self.major_dir) }
}

// --- CurveEval implementations ---

impl CurveEval for Line3 {
    fn point_at(&self, t: f64) -> DVec3 {
        self.origin + t * self.direction
    }
    fn tangent_at(&self, _t: f64) -> DVec3 {
        self.direction
    }
    fn derivative_at(&self, _t: f64) -> DVec3 {
        self.direction
    }
    fn derivative2_at(&self, _t: f64) -> DVec3 {
        DVec3::ZERO
    }
    fn derivative3_at(&self, _t: f64) -> DVec3 {
        DVec3::ZERO
    }
    fn curvature_at(&self, _t: f64) -> f64 {
        0.0
    }
    fn default_domain(&self) -> [f64; 2] {
        // OCCT Geom_Line::FirstParameter/LastParameter (Geom_Line.cxx L137/L144):
        // -Precision::Infinite() / Precision::Infinite().
        [-INFINITE_VALUE, INFINITE_VALUE]
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        -t
    }
}

impl CurveEval for Circle3 {
    fn point_at(&self, t: f64) -> DVec3 {
        self.center + self.x_dir * (self.radius * t.cos()) + self.y_dir * (self.radius * t.sin())
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        (-t.sin() * self.x_dir + t.cos() * self.y_dir).normalize()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        self.radius * (-t.sin() * self.x_dir + t.cos() * self.y_dir)
    }
    fn derivative2_at(&self, t: f64) -> DVec3 {
        // P''(t) = -R*(cos(t)·X + sin(t)·Y) = -(P(t) - center)
        -(self.point_at(t) - self.center)
    }
    fn derivative3_at(&self, t: f64) -> DVec3 {
        // P'''(t) = R*(sin(t)·X - cos(t)·Y) = -(1/R)*P'(t) = -derivative_at(t)
        -self.derivative_at(t)
    }
    fn curvature_at(&self, _t: f64) -> f64 {
        1.0 / self.radius
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
}

impl CurveEval for Ellipse3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let x_ax = self.major_dir;
        let y_ax = self.normal.cross(x_ax).normalize();
        self.center + self.major_radius * t.cos() * x_ax + self.minor_radius * t.sin() * y_ax
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let x_ax = self.major_dir;
        let y_ax = self.normal.cross(x_ax).normalize();
        (-self.major_radius * t.sin() * x_ax + self.minor_radius * t.cos() * y_ax).normalize()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        let x_ax = self.major_dir;
        let y_ax = self.normal.cross(x_ax).normalize();
        -self.major_radius * t.sin() * x_ax + self.minor_radius * t.cos() * y_ax
    }
    fn derivative2_at(&self, t: f64) -> DVec3 {
        let x_ax = self.major_dir;
        let y_ax = self.normal.cross(x_ax).normalize();
        -self.major_radius * t.cos() * x_ax - self.minor_radius * t.sin() * y_ax
    }
    fn derivative3_at(&self, t: f64) -> DVec3 {
        let x_ax = self.major_dir;
        let y_ax = self.normal.cross(x_ax).normalize();
        self.major_radius * t.sin() * x_ax - self.minor_radius * t.cos() * y_ax
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
}

impl CurveEval for Hyperbola3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let minor_dir = self.normal.cross(self.major_dir).normalize();
        self.center
            + self.semi_major * t.cosh() * self.major_dir
            + self.semi_minor * t.sinh() * minor_dir
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let minor_dir = self.normal.cross(self.major_dir).normalize();
        let v =
            self.semi_major * t.sinh() * self.major_dir + self.semi_minor * t.cosh() * minor_dir;
        v.normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        let minor_dir = self.normal.cross(self.major_dir).normalize();
        self.semi_major * t.sinh() * self.major_dir + self.semi_minor * t.cosh() * minor_dir
    }
    fn default_domain(&self) -> [f64; 2] {
        // OCCT Geom_Hyperbola::FirstParameter/LastParameter
        // (Geom_Hyperbola.cxx L94/L101): -Precision::Infinite() /
        // Precision::Infinite().
        [-INFINITE_VALUE, INFINITE_VALUE]
    }
}

impl CurveEval for Parabola3 {
    fn point_at(&self, t: f64) -> DVec3 {
        // OCCT Geom_Parabola (gp_Ax2 N, X): the cross-axis Y = N x X, so
        // dir_perp = normal x axis_dir forms the right-handed frame.
        let dir_perp = self.normal.cross(self.axis_dir).normalize();
        self.vertex + (t * t / (2.0 * self.focal_param)) * self.axis_dir + t * dir_perp
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let dir_perp = self.normal.cross(self.axis_dir).normalize();
        let v = (t / self.focal_param) * self.axis_dir + dir_perp;
        v.normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        let dir_perp = self.normal.cross(self.axis_dir).normalize();
        (t / self.focal_param) * self.axis_dir + dir_perp
    }
    fn default_domain(&self) -> [f64; 2] {
        // OCCT Geom_Parabola::FirstParameter/LastParameter
        // (Geom_Parabola.cxx L114/L128): -Precision::Infinite() /
        // Precision::Infinite().
        [-INFINITE_VALUE, INFINITE_VALUE]
    }
}

impl CurveEval for CircularHelix3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let axis = self.axis.normalize_or_zero();
        let mut x_axis = self.ref_dir - axis * self.ref_dir.dot(axis);
        if x_axis.length_squared() <= 1e-24 {
            x_axis = any_perpendicular(axis);
        } else {
            x_axis = x_axis.normalize();
        }
        let y_axis = axis.cross(x_axis).normalize_or_zero();
        let lead = self.pitch / (2.0 * PI);
        self.origin + self.radius * (t.cos() * x_axis + t.sin() * y_axis) + (lead * t) * axis
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let axis = self.axis.normalize_or_zero();
        let mut x_axis = self.ref_dir - axis * self.ref_dir.dot(axis);
        if x_axis.length_squared() <= 1e-24 {
            x_axis = any_perpendicular(axis);
        } else {
            x_axis = x_axis.normalize();
        }
        let y_axis = axis.cross(x_axis).normalize_or_zero();
        let lead = self.pitch / (2.0 * PI);
        (-self.radius * t.sin() * x_axis + self.radius * t.cos() * y_axis + lead * axis)
            .normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        let axis = self.axis.normalize_or_zero();
        let mut x_axis = self.ref_dir - axis * self.ref_dir.dot(axis);
        if x_axis.length_squared() <= 1e-24 {
            x_axis = any_perpendicular(axis);
        } else {
            x_axis = x_axis.normalize();
        }
        let y_axis = axis.cross(x_axis).normalize_or_zero();
        let lead = self.pitch / (2.0 * PI);
        -self.radius * t.sin() * x_axis + self.radius * t.cos() * y_axis + lead * axis
    }
    fn default_domain(&self) -> [f64; 2] {
        [-1e4, 1e4]
    }
}

impl CurveEval for SineWave3 {
    fn point_at(&self, t: f64) -> DVec3 {
        self.origin
            + t * self.baseline_dir
            + self.amplitude * (self.frequency * t + self.phase).sin() * self.amplitude_dir
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let v = self.baseline_dir
            + self.amplitude
                * self.frequency
                * (self.frequency * t + self.phase).cos()
                * self.amplitude_dir;
        v.normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        self.baseline_dir
            + self.amplitude
                * self.frequency
                * (self.frequency * t + self.phase).cos()
                * self.amplitude_dir
    }
    fn default_domain(&self) -> [f64; 2] {
        [-1e4, 1e4]
    }
}


impl Curve2dEval for Line2d {
    /// OCCT-aligned: P(t) = Location + t * Direction (Direction = gp_Dir2d = unit)
    fn point_at(&self, t: f64) -> DVec2 {
        self.origin + t * self.direction
    }
    /// OCCT-aligned: D1(t) = Direction = constant unit vector (gp_Dir2d invariant).
    fn tangent_at(&self, _t: f64) -> DVec2 {
        self.direction
    }
    fn derivative_at(&self, _t: f64) -> DVec2 {
        self.direction
    }
    /// OCCT-aligned: D2(t) = 0 for a line.
    fn derivative2_at(&self, _t: f64) -> DVec2 {
        DVec2::ZERO
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        -t
    }
    fn default_domain(&self) -> [f64; 2] {
        // OCCT Geom2d_Line::FirstParameter/LastParameter
        // (Geom2d_Line.cxx L144/L151): -Precision::Infinite() /
        // Precision::Infinite().
        [-INFINITE_VALUE, INFINITE_VALUE]
    }
}

impl Curve2dEval for Circle2d {
    fn point_at(&self, t: f64) -> DVec2 {
        // OCCT P(t) = Location + X_Dir * R*cos(t) + Y_Dir * R*sin(t)
        self.center + self.x_dir * (self.radius * t.cos()) + self.y_dir * (self.radius * t.sin())
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        (-t.sin() * self.x_dir + t.cos() * self.y_dir).normalize()
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        self.radius * (-t.sin() * self.x_dir + t.cos() * self.y_dir)
    }
    /// OCCT-aligned: D2(t) = -R * (cos(t)*X_Dir + sin(t)*Y_Dir).
    fn derivative2_at(&self, t: f64) -> DVec2 {
        self.radius * (-(t.cos()) * self.x_dir - t.sin() * self.y_dir)
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
}

impl Curve2dEval for Ellipse2d {
    fn point_at(&self, t: f64) -> DVec2 {
        self.center
            + self.major_dir * (self.major_radius * t.cos())
            + self.minor_dir * (self.minor_radius * t.sin())
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        (-self.major_radius * t.sin() * self.major_dir + self.minor_radius * t.cos() * self.minor_dir)
            .normalize()
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        -self.major_radius * t.sin() * self.major_dir + self.minor_radius * t.cos() * self.minor_dir
    }
    /// OCCT-aligned: D2(t) = -a·cos(t)·X_Dir - b·sin(t)·Y_Dir.
    fn derivative2_at(&self, t: f64) -> DVec2 {
        -self.major_radius * t.cos() * self.major_dir - self.minor_radius * t.sin() * self.minor_dir
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
    }
    fn is_closed(&self) -> bool {
        true
    }
    fn is_periodic(&self) -> bool {
        true
    }
    fn reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
}

impl Curve2dEval for Parabola2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let perp = DVec2::new(-self.axis_dir.y, self.axis_dir.x);
        self.origin + (t * t / (2.0 * self.focal_param)) * self.axis_dir + t * perp
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        let perp = DVec2::new(-self.axis_dir.y, self.axis_dir.x);
        (t / self.focal_param) * self.axis_dir + perp
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        // OCCT Geom2d_Parabola::FirstParameter/LastParameter
        // (Geom2d_Parabola.cxx L130/L137): -Precision::Infinite() /
        // Precision::Infinite().
        [-INFINITE_VALUE, INFINITE_VALUE]
    }
}

impl Curve2dEval for Hyperbola2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let minor = DVec2::new(-self.major_dir.y, self.major_dir.x);
        self.center
            + self.semi_major * t.cosh() * self.major_dir
            + self.semi_minor * t.sinh() * minor
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        let minor = DVec2::new(-self.major_dir.y, self.major_dir.x);
        self.semi_major * t.sinh() * self.major_dir + self.semi_minor * t.cosh() * minor
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        // OCCT Geom2d_Hyperbola::FirstParameter/LastParameter
        // (Geom2d_Hyperbola.cxx L147/L154): -Precision::Infinite() /
        // Precision::Infinite().
        [-INFINITE_VALUE, INFINITE_VALUE]
    }
}

impl Curve2dEval for CircleInvolute2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let r = self.base_radius.max(0.0);
        let x = r * (t.cos() + t * t.sin());
        let y = r * (t.sin() - t * t.cos());

        let ca = self.start_angle.cos();
        let sa = self.start_angle.sin();
        let xr = x * ca - y * sa;
        let yr = x * sa + y * ca;
        self.center + DVec2::new(xr, yr)
    }
}

impl Curve2dEval for ArchimedeanSpiral2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let r = self.a + self.b * t;
        let th = self.start_angle + t;
        self.center + DVec2::new(r * th.cos(), r * th.sin())
    }
}

impl Curve2dEval for LogarithmicSpiral2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let r = self.a * (self.b * t).exp();
        let th = self.start_angle + t;
        self.center + DVec2::new(r * th.cos(), r * th.sin())
    }
}

impl Curve2dEval for SineWave2d {
    fn point_at(&self, t: f64) -> DVec2 {
        DVec2::new(t, self.amplitude * (self.frequency * t + self.phase).sin())
    }
}
