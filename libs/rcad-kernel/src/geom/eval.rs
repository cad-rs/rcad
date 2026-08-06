use crate::geom::*;
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

// --- BoundedCurveEval implementations ---

impl BoundedCurveEval for BSplineCurve3 {
    fn degree(&self) -> usize { self.degree }
}

impl BoundedCurveEval for BezierCurve3 {
    fn degree(&self) -> usize {
        self.control_points.len().saturating_sub(1)
    }
}

// --- Curve3 type-group accessors (OCCT-aligned IsKind / DownCast equivalents) ---

impl Curve3 {
    /// Returns `true` if this curve is a conic (OCCT: `IsKind(Geom_Conic)`).
    pub fn is_conic(&self) -> bool {
        matches!(self, Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::Hyperbola(_) | Curve3::Parabola(_))
    }

    /// Returns `true` if this curve is bounded (OCCT: `IsKind(Geom_BoundedCurve)`).
    pub fn is_bounded(&self) -> bool {
        matches!(self, Curve3::BSpline(_) | Curve3::Bezier(_))
    }

    /// OCCT-aligned: downcast to conic trait object.
    pub fn as_conic(&self) -> Option<&dyn ConicEval> {
        match self {
            Curve3::Circle(c) => Some(c as &dyn ConicEval),
            Curve3::Ellipse(c) => Some(c as &dyn ConicEval),
            Curve3::Hyperbola(c) => Some(c as &dyn ConicEval),
            Curve3::Parabola(c) => Some(c as &dyn ConicEval),
            _ => None,
        }
    }

    /// OCCT-aligned: downcast to bounded curve trait object.
    pub fn as_bounded(&self) -> Option<&dyn BoundedCurveEval> {
        match self {
            Curve3::BSpline(c) => Some(c as &dyn BoundedCurveEval),
            Curve3::Bezier(c) => Some(c as &dyn BoundedCurveEval),
            _ => None,
        }
    }
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
    fn y_axis(&self) -> DVec2 { turn_2d(self.major_dir) }
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

// --- BoundedCurve2dEval implementations ---

impl BoundedCurve2dEval for BSplineCurve2 {
    fn degree(&self) -> usize { self.degree }
}

impl BoundedCurve2dEval for BezierCurve2 {
    fn degree(&self) -> usize {
        self.control_points.len().saturating_sub(1)
    }
}

// --- Curve2d type-group accessors (OCCT-aligned IsKind / DownCast equivalents) ---

impl Curve2d {
    /// Returns `true` if this curve is a conic (OCCT: `IsKind(Geom2d_Conic)`).
    pub fn is_conic(&self) -> bool {
        matches!(
            self,
            Curve2d::Circle(_) | Curve2d::Ellipse(_) | Curve2d::Hyperbola(_) | Curve2d::Parabola(_)
        )
    }

    /// Returns `true` if this curve is bounded (OCCT: `IsKind(Geom2d_BoundedCurve)`).
    pub fn is_bounded(&self) -> bool {
        matches!(self, Curve2d::BSpline(_) | Curve2d::Bezier(_))
    }

    /// OCCT-aligned: downcast to 2D conic trait object.
    pub fn as_conic(&self) -> Option<&dyn Conic2dEval> {
        match self {
            Curve2d::Circle(c) => Some(c as &dyn Conic2dEval),
            Curve2d::Ellipse(c) => Some(c as &dyn Conic2dEval),
            Curve2d::Hyperbola(c) => Some(c as &dyn Conic2dEval),
            Curve2d::Parabola(c) => Some(c as &dyn Conic2dEval),
            _ => None,
        }
    }

    /// OCCT-aligned: downcast to bounded 2D curve trait object.
    pub fn as_bounded(&self) -> Option<&dyn BoundedCurve2dEval> {
        match self {
            Curve2d::BSpline(c) => Some(c as &dyn BoundedCurve2dEval),
            Curve2d::Bezier(c) => Some(c as &dyn BoundedCurve2dEval),
            _ => None,
        }
    }
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
        [f64::NEG_INFINITY, f64::INFINITY]
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
        [-1e4, 1e4] // unbounded; caller trims as needed
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
        [-1e4, 1e4] // unbounded
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

impl CurveEval for TrimmedCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        self.curve.point_at(self.map_param(t))
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        self.curve.tangent_at(self.map_param(t))
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        self.curve.derivative_at(self.map_param(t))
    }
    fn default_domain(&self) -> [f64; 2] {
        [self.first, self.last]
    }
}

impl CurveEval for Curve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.point_at(t),
            Curve3::Circle(c) => c.point_at(t),
            Curve3::Ellipse(c) => c.point_at(t),
            Curve3::BSpline(c) => c.point_at(t),
            Curve3::Bezier(c) => c.point_at(t),
            Curve3::Offset(c) => c.point_at(t),
            Curve3::Hyperbola(c) => c.point_at(t),
            Curve3::Parabola(c) => c.point_at(t),
            Curve3::CircularHelix(c) => c.point_at(t),
            Curve3::SineWave(c) => c.point_at(t),
            Curve3::Trimmed(tc) => tc.point_at(t),
        }
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.tangent_at(t),
            Curve3::Circle(c) => c.tangent_at(t),
            Curve3::Ellipse(c) => c.tangent_at(t),
            Curve3::BSpline(c) => c.tangent_at(t),
            Curve3::Bezier(c) => c.tangent_at(t),
            Curve3::Offset(c) => c.tangent_at(t),
            Curve3::Hyperbola(c) => c.tangent_at(t),
            Curve3::Parabola(c) => c.tangent_at(t),
            Curve3::CircularHelix(c) => c.tangent_at(t),
            Curve3::SineWave(c) => c.tangent_at(t),
            Curve3::Trimmed(tc) => tc.tangent_at(t),
        }
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.derivative_at(t),
            Curve3::Circle(c) => c.derivative_at(t),
            Curve3::Ellipse(c) => c.derivative_at(t),
            Curve3::BSpline(c) => c.derivative_at(t),
            Curve3::Bezier(c) => c.derivative_at(t),
            Curve3::Offset(c) => c.derivative_at(t),
            Curve3::Hyperbola(c) => c.derivative_at(t),
            Curve3::Parabola(c) => c.derivative_at(t),
            Curve3::CircularHelix(c) => c.derivative_at(t),
            Curve3::SineWave(c) => c.derivative_at(t),
            Curve3::Trimmed(tc) => tc.derivative_at(t),
        }
    }
    fn derivative2_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.derivative2_at(t),
            Curve3::Circle(c) => c.derivative2_at(t),
            Curve3::Ellipse(c) => c.derivative2_at(t),
            Curve3::BSpline(c) => c.derivative2_at(t),
            Curve3::Bezier(c) => c.derivative2_at(t),
            Curve3::Offset(c) => c.derivative2_at(t),
            Curve3::Hyperbola(c) => c.derivative2_at(t),
            Curve3::Parabola(c) => c.derivative2_at(t),
            Curve3::CircularHelix(c) => c.derivative2_at(t),
            Curve3::SineWave(c) => c.derivative2_at(t),
            Curve3::Trimmed(tc) => tc.derivative2_at(t),
        }
    }
    fn derivative3_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.derivative3_at(t),
            Curve3::Circle(c) => c.derivative3_at(t),
            Curve3::Ellipse(c) => c.derivative3_at(t),
            Curve3::BSpline(c) => c.derivative3_at(t),
            Curve3::Bezier(c) => c.derivative3_at(t),
            Curve3::Offset(c) => c.derivative3_at(t),
            Curve3::Hyperbola(c) => c.derivative3_at(t),
            Curve3::Parabola(c) => c.derivative3_at(t),
            Curve3::CircularHelix(c) => c.derivative3_at(t),
            Curve3::SineWave(c) => c.derivative3_at(t),
            Curve3::Trimmed(tc) => tc.derivative3_at(t),
        }
    }
    fn curvature_at(&self, t: f64) -> f64 {
        match self {
            Curve3::Line(c) => c.curvature_at(t),
            Curve3::Circle(c) => c.curvature_at(t),
            Curve3::Ellipse(c) => c.curvature_at(t),
            Curve3::BSpline(c) => c.curvature_at(t),
            Curve3::Bezier(c) => c.curvature_at(t),
            Curve3::Offset(c) => c.curvature_at(t),
            Curve3::Hyperbola(c) => c.curvature_at(t),
            Curve3::Parabola(c) => c.curvature_at(t),
            Curve3::CircularHelix(c) => c.curvature_at(t),
            Curve3::SineWave(c) => c.curvature_at(t),
            Curve3::Trimmed(tc) => tc.curvature_at(t),
        }
    }
    fn transformed_parameter(&self, t: f64) -> f64 {
        match self {
            Curve3::Line(c) => c.transformed_parameter(t),
            Curve3::Circle(c) => c.transformed_parameter(t),
            Curve3::Ellipse(c) => c.transformed_parameter(t),
            Curve3::BSpline(c) => c.transformed_parameter(t),
            Curve3::Bezier(c) => c.transformed_parameter(t),
            Curve3::Offset(c) => c.transformed_parameter(t),
            Curve3::Hyperbola(c) => c.transformed_parameter(t),
            Curve3::Parabola(c) => c.transformed_parameter(t),
            Curve3::CircularHelix(c) => c.transformed_parameter(t),
            Curve3::SineWave(c) => c.transformed_parameter(t),
            Curve3::Trimmed(tc) => tc.transformed_parameter(t),
        }
    }
    fn parametric_transformation(&self) -> f64 {
        match self {
            Curve3::Line(c) => c.parametric_transformation(),
            Curve3::Circle(c) => c.parametric_transformation(),
            Curve3::Ellipse(c) => c.parametric_transformation(),
            Curve3::BSpline(c) => c.parametric_transformation(),
            Curve3::Bezier(c) => c.parametric_transformation(),
            Curve3::Offset(c) => c.parametric_transformation(),
            Curve3::Hyperbola(c) => c.parametric_transformation(),
            Curve3::Parabola(c) => c.parametric_transformation(),
            Curve3::CircularHelix(c) => c.parametric_transformation(),
            Curve3::SineWave(c) => c.parametric_transformation(),
            Curve3::Trimmed(tc) => tc.parametric_transformation(),
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        match self {
            Curve3::Line(c) => c.default_domain(),
            Curve3::Circle(c) => c.default_domain(),
            Curve3::Ellipse(c) => c.default_domain(),
            Curve3::BSpline(c) => c.default_domain(),
            Curve3::Bezier(c) => c.default_domain(),
            Curve3::Offset(c) => c.default_domain(),
            Curve3::Hyperbola(c) => c.default_domain(),
            Curve3::Parabola(c) => c.default_domain(),
            Curve3::CircularHelix(c) => c.default_domain(),
            Curve3::SineWave(c) => c.default_domain(),
            Curve3::Trimmed(tc) => tc.default_domain(),
        }
    }
}

// --- SurfaceEval implementations ---

impl SurfaceEval for Plane {
    /// OCCT-aligned: P(u,v) = origin + u*u_dir + v*v_dir using stored axes.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.origin + u * self.u_dir + v * self.v_dir
    }
    fn normal_at(&self, _u: f64, _v: f64) -> DVec3 {
        self.normal
    }
    fn default_domain(&self) -> [f64; 4] {
        [
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
        ]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        (
            self.origin + u * self.u_dir + v * self.v_dir,
            self.u_dir,
            self.v_dir,
        )
    }
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        (
            self.origin + u * self.u_dir + v * self.v_dir,
            self.u_dir,
            self.v_dir,
            DVec3::ZERO, // Puu
            DVec3::ZERO, // Puv
            DVec3::ZERO, // Pvv
        )
    }
}

impl SurfaceEval for CylindricalSurface {
    /// u = azimuth angle [0, 2π], v = height along axis.
    /// OCCT-aligned: uses stored ref_dir for deterministic UV mapping.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.axis.cross(x_ax).normalize();
        self.origin + self.radius * (u.cos() * x_ax + u.sin() * y_ax) + v * self.axis
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.axis.cross(x_ax).normalize();
        (u.cos() * x_ax + u.sin() * y_ax).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, f64::NEG_INFINITY, f64::INFINITY]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let p = self.origin + self.radius * (cu * x_ax + su * y_ax) + v * self.axis;
        let dpu = self.radius * (-su * x_ax + cu * y_ax);
        (p, dpu, self.axis)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
}

impl CylindricalSurface {
    /// UV coordinates of world point `p` relative to this cylindrical surface.
    ///
    /// `u` = azimuth (−π, π], `v` = height along axis, matching
    /// [`SurfaceEval::point_at`] (which uses the stored `ref_dir` for u=0).
    /// Off-surface points are radially projected onto the cylinder.
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        let axis = self.axis.normalize_or_zero();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let d = p - self.origin;
        let v = d.dot(axis);
        let radial = d - axis * v;
        let r = radial.length();
        if r < 1e-15 {
            return DVec2::new(0.0, v);
        }
        let u = radial.dot(y_ax).atan2(radial.dot(x_ax));
        DVec2::new(u, v)
    }
}

impl SphericalSurface {
    /// Construct a sphere with `ref_dir` derived from [`any_perpendicular(axis)`](any_perpendicular).
    pub fn new(center: Point3, axis: Vec3, radius: f64) -> Self {
        Self {
            center,
            axis,
            radius,
            ref_dir: any_perpendicular(axis),
        }
    }

    /// Construct a sphere with an explicit `ref_dir` (used after mirroring / transforming).
    pub fn new_with_ref_dir(center: Point3, axis: Vec3, radius: f64, ref_dir: Vec3) -> Self {
        Self {
            center,
            axis,
            radius,
            ref_dir,
        }
    }

    /// Spherical coordinates of world point `p`: longitude `u` ∈ (−π, π], colatitude `v` ∈ [0, π],
    /// matching [`SurfaceEval::point_at`] / `properties` sphere helpers (radial projection when `p`
    /// is off the surface).
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        let ax = self.axis.normalize_or_zero();
        let r = self.radius;
        if r < 1e-15 {
            return DVec2::ZERO;
        }
        let w = (p - self.center) / r;
        if w.length_squared() < 1e-20 {
            return DVec2::ZERO;
        }
        let w = w.normalize();
        let v = w.dot(ax).clamp(-1.0, 1.0).acos();
        let x_ax = self.ref_dir.normalize();
        let y_ax = ax.cross(x_ax).normalize();
        let w_t = w - ax * w.dot(ax);
        if w_t.length_squared() < 1e-12 {
            return DVec2::new(0.0, v);
        }
        let w_t = w_t.normalize();
        let u = w_t.dot(y_ax).atan2(w_t.dot(x_ax));
        DVec2::new(u, v)
    }
}

impl SurfaceEval for SphericalSurface {
    /// u = longitude [0, 2π], v = colatitude [0, π] (0 = north pole).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize();
        let y_ax = self.axis.cross(x_ax).normalize();
        self.center
            + self.radius * (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * self.axis)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let p = self.point_at(u, v);
        (p - self.center).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, PI]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let x_ax = self.ref_dir.normalize();
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let radial = cu * x_ax + su * y_ax;
        let p = self.center + self.radius * (sv * radial + cv * self.axis);
        let dpu = self.radius * sv * (-su * x_ax + cu * y_ax);
        let dpv = self.radius * (cv * radial - sv * self.axis);
        (p, dpu, dpv)
    }
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        // OCCT Geom_SphericalSurface::D2. In the basis (X, Y, A) with
        // radial = cos(u)X + sin(u)Y:
        //   Puu = R sin(v)(-cos(u)X - sin(u)Y)
        //   Puv = R cos(v)(-sin(u)X + cos(u)Y)
        //   Pvv = -(P - center)
        let x_ax = self.ref_dir.normalize();
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let radial = cu * x_ax + su * y_ax;
        let p = self.center + self.radius * (sv * radial + cv * self.axis);
        let dpu = self.radius * sv * (-su * x_ax + cu * y_ax);
        let dpv = self.radius * (cv * radial - sv * self.axis);
        let d2u = self.radius * sv * (-cu * x_ax - su * y_ax);
        let duv = self.radius * cv * (-su * x_ax + cu * y_ax);
        let d2v = -(p - self.center);
        (p, dpu, dpv, d2u, duv, d2v)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
    fn v_reversed_parameter(&self, t: f64) -> f64 {
        PI - t
    } // OCCT: colatitude [0, π]
}

impl SurfaceEval for ConicalSurface {
    /// u = azimuth [0, 2π], v = distance along the cone generatrix from the
    /// reference circle at `self.apex`.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let radial = self.radius_at_slant(v);
        let axial = self.axial_from_slant(v);
        self.apex + axial * axis + radial * (u.cos() * x_ax + u.sin() * y_ax)
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let radial = u.cos() * x_ax + u.sin() * y_ax;
        let half = self.half_angle_rad;
        (radial * half.cos() - axis * half.sin()).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, f64::INFINITY]
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let radial = self.radius_at_slant(v);
        let axial = self.axial_from_slant(v);
        // d(radius)/dv = sin(half_angle), d(axial)/dv = cos(half_angle)
        let half = self.half_angle_rad;
        let dr = half.sin();
        let da = half.cos();
        let r_vec = cu * x_ax + su * y_ax;
        let p = self.apex + axial * axis + radial * r_vec;
        let dpu = radial * (-su * x_ax + cu * y_ax);
        let dpv = da * axis + dr * r_vec;
        (p, dpu, dpv)
    }
}

impl SurfaceEval for ToroidalSurface {
    /// u = major angle [0, 2π], v = minor angle [0, 2π].
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        let tube_center = self.center + self.major_radius * (u.cos() * x_ax + u.sin() * y_ax);
        let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
        tube_center + self.minor_radius * (v.cos() * radial + v.sin() * self.axis)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
        (v.cos() * radial + v.sin() * self.axis).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, 2.0 * PI]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let r_vec = cu * x_ax + su * y_ax;
        let r_perp = -su * x_ax + cu * y_ax;
        let r_major = self.major_radius;
        let r_minor = self.minor_radius;
        let tube = r_major + r_minor * cv;
        let p = self.center + tube * r_vec + r_minor * sv * self.axis;
        let dpu = tube * r_perp;
        let dpv = -r_minor * sv * r_vec + r_minor * cv * self.axis;
        (p, dpu, dpv)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
    fn is_v_closed(&self) -> bool {
        true
    }
    fn is_v_periodic(&self) -> bool {
        true
    }
    fn v_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
}

impl ToroidalSurface {
    /// UV coordinates of world point `p` relative to this toroidal surface.
    ///
    /// `u` = major angle (−π, π], `v` = minor angle [0, 2π),
    /// matching [`SurfaceEval::point_at`].  When `p` is on the surface
    /// the returned `(u, v)` is exact; off-surface points project onto
    /// the tube center circle in the radial direction.
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        use std::f64::consts::TAU;
        let axis = self.axis.normalize_or_zero();
        let x_ax = any_perpendicular(axis);
        let y_ax = axis.cross(x_ax).normalize();
        let local = p - self.center;
        let axial = local.dot(axis);
        let radial_vec = local - axis * axial;
        let radial_dist = radial_vec.length();

        // u = azimuth around main axis
        let u = if radial_dist < 1e-15 {
            0.0
        } else {
            let rn = radial_vec / radial_dist;
            rn.dot(y_ax).atan2(rn.dot(x_ax))
        };

        // v = angle around tube:
        // On surface: radial_dist = R + r·cos(v), axial = r·sin(v)
        //   → v = atan2(axial, radial_dist - R)
        let v_base = axial.atan2(radial_dist - self.major_radius);
        // Convert v from [-π, π] to [0, 2π)
        let v = if v_base < 0.0 { v_base + TAU } else { v_base };

        DVec2::new(u, v)
    }
}

impl SurfaceEval for EllipsoidalSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        self.center
            + self.radius_x * v.sin() * u.cos() * x_axis
            + self.radius_y * v.sin() * u.sin() * y_axis
            + self.radius_z * v.cos() * axis
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let p = self.point_at(u, v) - self.center;
        let x = p.dot(x_axis);
        let y = p.dot(y_axis);
        let z = p.dot(axis);
        let grad = (x / (self.radius_x * self.radius_x)) * x_axis
            + (y / (self.radius_y * self.radius_y)) * y_axis
            + (z / (self.radius_z * self.radius_z)) * axis;
        grad.normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, PI]
    }
}

impl SurfaceEval for HelicoidSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let lead = self.pitch / (2.0 * PI);
        self.origin + v * (u.cos() * x_axis + u.sin() * y_axis) + (lead * u) * axis
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let lead = self.pitch / (2.0 * PI);
        let du = v * (-u.sin() * x_axis + u.cos() * y_axis) + lead * axis;
        let dv = u.cos() * x_axis + u.sin() * y_axis;
        du.cross(dv).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [-2.0 * PI, 2.0 * PI, -10.0, 10.0]
    }
}

impl SurfaceEval for LinearExtrusionSurface {
    /// u = profile parameter, v = extrusion distance along direction.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.profile.point_at(u) + v * self.direction
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let tangent = self.profile.tangent_at(u);
        let n = tangent.cross(self.direction);
        if n.length_squared() < 1e-20 {
            return DVec3::Z;
        }
        n.normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        let [t1, t2] = self.profile.default_domain();
        [t1, t2, f64::NEG_INFINITY, f64::INFINITY]
    }
}

impl SurfaceEval for RevolutionSurface {
    /// u = azimuth angle [0, 2π], v = profile parameter.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let p = self.profile.point_at(v);
        let d = p - self.axis_origin;
        let d_par = self.axis_dir * d.dot(self.axis_dir);
        let d_perp = d - d_par;
        self.axis_origin + d_par + d_perp * u.cos() + self.axis_dir.cross(d_perp) * u.sin()
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        let n = du.cross(dv);
        if n.length_squared() < 1e-20 {
            return DVec3::Z;
        }
        n.normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        let [t1, t2] = self.profile.default_domain();
        [0.0, 2.0 * PI, t1, t2]
    }
}

impl SurfaceEval for TrimmedSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.basis.point_at(u, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        self.basis.normal_at(u, v)
    }
    fn default_domain(&self) -> [f64; 4] {
        self.trim
    }
}

impl SurfaceEval for Surface3 {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface3::Plane(s) => s.point_at(u, v),
            Surface3::Cylinder(s) => s.point_at(u, v),
            Surface3::Sphere(s) => s.point_at(u, v),
            Surface3::Cone(s) => s.point_at(u, v),
            Surface3::Torus(s) => s.point_at(u, v),
            Surface3::Ellipsoid(s) => s.point_at(u, v),
            Surface3::Helicoid(s) => s.point_at(u, v),
            Surface3::Pipe(s) => s.point_at(u, v),
            Surface3::BSpline(s) => s.point_at(u, v),
            Surface3::LinearExtrusion(s) => s.point_at(u, v),
            Surface3::Revolution(s) => s.point_at(u, v),
            Surface3::Ruled(s) => s.point_at(u, v),
            Surface3::Coons(s) => s.point_at(u, v),
            Surface3::Bezier(s) => s.point_at(u, v),
            Surface3::TriBezier(s) => s.point_at(u, v),
            Surface3::Offset(s) => s.point_at(u, v),
            Surface3::Trimmed(s) => s.point_at(u, v),
        }
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface3::Plane(s) => s.normal_at(u, v),
            Surface3::Cylinder(s) => s.normal_at(u, v),
            Surface3::Sphere(s) => s.normal_at(u, v),
            Surface3::Cone(s) => s.normal_at(u, v),
            Surface3::Torus(s) => s.normal_at(u, v),
            Surface3::Ellipsoid(s) => s.normal_at(u, v),
            Surface3::Helicoid(s) => s.normal_at(u, v),
            Surface3::Pipe(s) => s.normal_at(u, v),
            Surface3::BSpline(s) => s.normal_at(u, v),
            Surface3::LinearExtrusion(s) => s.normal_at(u, v),
            Surface3::Revolution(s) => s.normal_at(u, v),
            Surface3::Ruled(s) => s.normal_at(u, v),
            Surface3::Coons(s) => s.normal_at(u, v),
            Surface3::Bezier(s) => s.normal_at(u, v),
            Surface3::TriBezier(s) => s.normal_at(u, v),
            Surface3::Offset(s) => s.normal_at(u, v),
            Surface3::Trimmed(s) => s.normal_at(u, v),
        }
    }
    fn default_domain(&self) -> [f64; 4] {
        match self {
            Surface3::Plane(s) => s.default_domain(),
            Surface3::Cylinder(s) => s.default_domain(),
            Surface3::Sphere(s) => s.default_domain(),
            Surface3::Cone(s) => s.default_domain(),
            Surface3::Torus(s) => s.default_domain(),
            Surface3::Ellipsoid(s) => s.default_domain(),
            Surface3::Helicoid(s) => s.default_domain(),
            Surface3::Pipe(s) => s.default_domain(),
            Surface3::BSpline(s) => s.default_domain(),
            Surface3::LinearExtrusion(s) => s.default_domain(),
            Surface3::Revolution(s) => s.default_domain(),
            Surface3::Ruled(s) => s.default_domain(),
            Surface3::Coons(s) => s.default_domain(),
            Surface3::Bezier(s) => s.default_domain(),
            Surface3::TriBezier(s) => s.default_domain(),
            Surface3::Offset(s) => s.default_domain(),
            Surface3::Trimmed(s) => s.default_domain(),
        }
    }
    fn is_u_closed(&self) -> bool {
        // OCCT Geom_Surface::IsUClosed — elementary surfaces of revolution are
        // closed in U (cylinder/cone/sphere/torus), planes and others are not.
        match self {
            Surface3::Plane(_) => false,
            Surface3::Cylinder(s) => s.is_u_closed(),
            Surface3::Sphere(s) => s.is_u_closed(),
            Surface3::Cone(s) => s.is_u_closed(),
            Surface3::Torus(s) => s.is_u_closed(),
            Surface3::Ellipsoid(_) => false,
            Surface3::Helicoid(_) => false,
            Surface3::Pipe(_) => false,
            Surface3::BSpline(_) => false,
            Surface3::LinearExtrusion(_) => false,
            Surface3::Revolution(_) => false,
            Surface3::Ruled(_) => false,
            Surface3::Coons(_) => false,
            Surface3::Bezier(_) => false,
            Surface3::TriBezier(_) => false,
            Surface3::Offset(_) => false,
            Surface3::Trimmed(_) => false,
        }
    }
    fn is_v_closed(&self) -> bool {
        // OCCT Geom_Surface::IsVClosed — only the torus is closed in V.
        match self {
            Surface3::Plane(_) => false,
            Surface3::Cylinder(s) => s.is_v_closed(),
            Surface3::Sphere(s) => s.is_v_closed(),
            Surface3::Cone(s) => s.is_v_closed(),
            Surface3::Torus(s) => s.is_v_closed(),
            Surface3::Ellipsoid(_) => false,
            Surface3::Helicoid(_) => false,
            Surface3::Pipe(_) => false,
            Surface3::BSpline(_) => false,
            Surface3::LinearExtrusion(_) => false,
            Surface3::Revolution(_) => false,
            Surface3::Ruled(_) => false,
            Surface3::Coons(_) => false,
            Surface3::Bezier(_) => false,
            Surface3::TriBezier(_) => false,
            Surface3::Offset(_) => false,
            Surface3::Trimmed(_) => false,
        }
    }
    fn is_u_periodic(&self) -> bool {
        match self {
            Surface3::Plane(s) => s.is_u_periodic(),
            Surface3::Cylinder(s) => s.is_u_periodic(),
            Surface3::Sphere(s) => s.is_u_periodic(),
            Surface3::Cone(s) => s.is_u_periodic(),
            Surface3::Torus(s) => s.is_u_periodic(),
            Surface3::Ellipsoid(s) => s.is_u_periodic(),
            Surface3::Helicoid(s) => s.is_u_periodic(),
            Surface3::Pipe(s) => s.is_u_periodic(),
            Surface3::BSpline(s) => s.is_u_periodic(),
            Surface3::LinearExtrusion(s) => s.is_u_periodic(),
            Surface3::Revolution(s) => s.is_u_periodic(),
            Surface3::Ruled(s) => s.is_u_periodic(),
            Surface3::Coons(s) => s.is_u_periodic(),
            Surface3::Bezier(s) => s.is_u_periodic(),
            Surface3::TriBezier(s) => s.is_u_periodic(),
            Surface3::Offset(s) => s.is_u_periodic(),
            Surface3::Trimmed(s) => s.is_u_periodic(),
        }
    }
    fn is_v_periodic(&self) -> bool {
        match self {
            Surface3::Plane(s) => s.is_v_periodic(),
            Surface3::Cylinder(s) => s.is_v_periodic(),
            Surface3::Sphere(s) => s.is_v_periodic(),
            Surface3::Cone(s) => s.is_v_periodic(),
            Surface3::Torus(s) => s.is_v_periodic(),
            Surface3::Ellipsoid(s) => s.is_v_periodic(),
            Surface3::Helicoid(s) => s.is_v_periodic(),
            Surface3::Pipe(s) => s.is_v_periodic(),
            Surface3::BSpline(s) => s.is_v_periodic(),
            Surface3::LinearExtrusion(s) => s.is_v_periodic(),
            Surface3::Revolution(s) => s.is_v_periodic(),
            Surface3::Ruled(s) => s.is_v_periodic(),
            Surface3::Coons(s) => s.is_v_periodic(),
            Surface3::Bezier(s) => s.is_v_periodic(),
            Surface3::TriBezier(s) => s.is_v_periodic(),
            Surface3::Offset(s) => s.is_v_periodic(),
            Surface3::Trimmed(s) => s.is_v_periodic(),
        }
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        match self {
            Surface3::Plane(s) => s.derivatives(u, v),
            Surface3::Cylinder(s) => s.derivatives(u, v),
            Surface3::Sphere(s) => s.derivatives(u, v),
            Surface3::Cone(s) => s.derivatives(u, v),
            Surface3::Torus(s) => s.derivatives(u, v),
            Surface3::Ellipsoid(s) => s.derivatives(u, v),
            Surface3::Helicoid(s) => s.derivatives(u, v),
            Surface3::Pipe(s) => s.derivatives(u, v),
            Surface3::BSpline(s) => s.derivatives(u, v),
            Surface3::LinearExtrusion(s) => s.derivatives(u, v),
            Surface3::Revolution(s) => s.derivatives(u, v),
            Surface3::Ruled(s) => s.derivatives(u, v),
            Surface3::Coons(s) => s.derivatives(u, v),
            Surface3::Bezier(s) => s.derivatives(u, v),
            Surface3::TriBezier(s) => s.derivatives(u, v),
            Surface3::Offset(s) => s.derivatives(u, v),
            Surface3::Trimmed(s) => s.derivatives(u, v),
        }
    }
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        match self {
            Surface3::Plane(s) => s.derivatives2(u, v),
            Surface3::Cylinder(s) => s.derivatives2(u, v),
            Surface3::Sphere(s) => s.derivatives2(u, v),
            Surface3::Cone(s) => s.derivatives2(u, v),
            Surface3::Torus(s) => s.derivatives2(u, v),
            Surface3::Ellipsoid(s) => s.derivatives2(u, v),
            Surface3::Helicoid(s) => s.derivatives2(u, v),
            Surface3::Pipe(s) => s.derivatives2(u, v),
            Surface3::BSpline(s) => s.derivatives2(u, v),
            Surface3::LinearExtrusion(s) => s.derivatives2(u, v),
            Surface3::Revolution(s) => s.derivatives2(u, v),
            Surface3::Ruled(s) => s.derivatives2(u, v),
            Surface3::Coons(s) => s.derivatives2(u, v),
            Surface3::Bezier(s) => s.derivatives2(u, v),
            Surface3::TriBezier(s) => s.derivatives2(u, v),
            Surface3::Offset(s) => s.derivatives2(u, v),
            Surface3::Trimmed(s) => s.derivatives2(u, v),
        }
    }
}

/// OCCT-aligned: `Geom_ElementarySurface` intermediate abstract class.
///
// --- BoundedSurfaceEval implementations ---

impl BoundedSurfaceEval for BSplineSurface {
    fn degree_u(&self) -> usize { self.degree_u }
    fn degree_v(&self) -> usize { self.degree_v }
}

impl BoundedSurfaceEval for BezierSurface {
    fn degree_u(&self) -> usize {
        self.control_points.len().saturating_sub(1)
    }
    fn degree_v(&self) -> usize {
        self.control_points.first().map_or(0, |r| r.len().saturating_sub(1))
    }
}

// --- SweptSurfaceEval implementations ---

impl SweptSurfaceEval for LinearExtrusionSurface {
    fn profile(&self) -> &Curve3 { &self.profile }
}

impl SweptSurfaceEval for RevolutionSurface {
    fn profile(&self) -> &Curve3 { &self.profile }
}

// --- ElementarySurfaceEval implementations ---

impl ElementarySurfaceEval for Plane {
    fn position(&self) -> DVec3 { self.origin }
    fn axis_dir(&self) -> DVec3 { self.normal }
    fn x_axis(&self) -> DVec3 { self.u_dir }
    fn y_axis(&self) -> DVec3 { self.v_dir }
}

impl ElementarySurfaceEval for CylindricalSurface {
    fn position(&self) -> DVec3 { self.origin }
    fn axis_dir(&self) -> DVec3 { self.axis }
    fn x_axis(&self) -> DVec3 { self.ref_dir.normalize_or_zero() }
    fn y_axis(&self) -> DVec3 { self.axis.cross(self.ref_dir).normalize_or_zero() }
}

impl ElementarySurfaceEval for SphericalSurface {
    fn position(&self) -> DVec3 { self.center }
    fn axis_dir(&self) -> DVec3 { self.axis }
    fn x_axis(&self) -> DVec3 { self.ref_dir.normalize_or_zero() }
    fn y_axis(&self) -> DVec3 { self.axis.cross(self.ref_dir).normalize_or_zero() }
}

impl ElementarySurfaceEval for ConicalSurface {
    fn position(&self) -> DVec3 { self.apex_point() }
    fn axis_dir(&self) -> DVec3 { self.axis_dir() }
    fn x_axis(&self) -> DVec3 { self.ref_dir.normalize_or_zero() }
    fn y_axis(&self) -> DVec3 { self.axis_dir().cross(self.ref_dir.normalize_or_zero()).normalize_or_zero() }
}

impl ElementarySurfaceEval for ToroidalSurface {
    fn position(&self) -> DVec3 { self.center }
    fn axis_dir(&self) -> DVec3 { self.axis }
    fn x_axis(&self) -> DVec3 { any_perpendicular(self.axis) }
    fn y_axis(&self) -> DVec3 { self.axis.cross(any_perpendicular(self.axis)).normalize_or_zero() }
}

// --- Surface3 type-group accessors ---

impl Surface3 {
    /// Returns `true` if this surface is an elementary surface
    /// (OCCT: `IsKind(Geom_ElementarySurface)`).
    pub fn is_elementary(&self) -> bool {
        matches!(
            self,
            Surface3::Plane(_)
                | Surface3::Cylinder(_)
                | Surface3::Sphere(_)
                | Surface3::Cone(_)
                | Surface3::Torus(_)
        )
    }

    /// Returns `true` if this surface is bounded
    /// (OCCT: `IsKind(Geom_BoundedSurface)`).
    pub fn is_bounded(&self) -> bool {
        matches!(self, Surface3::BSpline(_) | Surface3::Bezier(_))
    }

    /// OCCT-aligned: downcast to elementary surface trait object.
    pub fn as_elementary(&self) -> Option<&dyn ElementarySurfaceEval> {
        match self {
            Surface3::Plane(s) => Some(s as &dyn ElementarySurfaceEval),
            Surface3::Cylinder(s) => Some(s as &dyn ElementarySurfaceEval),
            Surface3::Sphere(s) => Some(s as &dyn ElementarySurfaceEval),
            Surface3::Cone(s) => Some(s as &dyn ElementarySurfaceEval),
            Surface3::Torus(s) => Some(s as &dyn ElementarySurfaceEval),
            _ => None,
        }
    }

    /// OCCT-aligned: downcast to bounded surface trait object.
    pub fn as_bounded(&self) -> Option<&dyn BoundedSurfaceEval> {
        match self {
            Surface3::BSpline(s) => Some(s as &dyn BoundedSurfaceEval),
            Surface3::Bezier(s) => Some(s as &dyn BoundedSurfaceEval),
            _ => None,
        }
    }

    /// Returns `true` if this surface is a swept surface
    /// (OCCT: `IsKind(Geom_SweptSurface)`).
    pub fn is_swept(&self) -> bool {
        matches!(self, Surface3::LinearExtrusion(_) | Surface3::Revolution(_))
    }

    /// OCCT-aligned: downcast to swept surface trait object.
    pub fn as_swept(&self) -> Option<&dyn SweptSurfaceEval> {
        match self {
            Surface3::LinearExtrusion(s) => Some(s as &dyn SweptSurfaceEval),
            Surface3::Revolution(s) => Some(s as &dyn SweptSurfaceEval),
            _ => None,
        }
    }
}

fn remap_unit_to_curve_domain(curve: &Curve3, t: f64) -> f64 {
    let [t0, t1] = curve.default_domain();
    if !t0.is_finite() || !t1.is_finite() {
        return t;
    }
    t0 + (t1 - t0) * t
}

fn projected_frame_from_tangent(tangent: DVec3, ref_dir: DVec3) -> (DVec3, DVec3) {
    let tangent = tangent.normalize_or_zero();
    let mut x_axis = ref_dir - tangent * ref_dir.dot(tangent);
    if x_axis.length_squared() <= 1e-24 {
        x_axis = any_perpendicular(tangent);
    } else {
        x_axis = x_axis.normalize();
    }
    let y_axis = tangent.cross(x_axis).normalize_or_zero();
    (x_axis, y_axis)
}

impl SurfaceEval for PipeSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let center = self.spine.point_at(v);
        let tangent = self.spine.tangent_at(v);
        let (x_axis, y_axis) = projected_frame_from_tangent(tangent, self.ref_dir);
        center + self.radius * (u.cos() * x_axis + u.sin() * y_axis)
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = self.point_at(u + eps, v) - self.point_at(u - eps, v);
        let dv = self.point_at(u, v + eps) - self.point_at(u, v - eps);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        let [v0, v1] = self.spine.default_domain();
        [0.0, 2.0 * PI, v0, v1]
    }
}

impl SurfaceEval for CoonsSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let south = self
            .south
            .point_at(remap_unit_to_curve_domain(&self.south, u));
        let north = self
            .north
            .point_at(remap_unit_to_curve_domain(&self.north, u));
        let west = self
            .west
            .point_at(remap_unit_to_curve_domain(&self.west, v));
        let east = self
            .east
            .point_at(remap_unit_to_curve_domain(&self.east, v));

        let p00 = self
            .south
            .point_at(remap_unit_to_curve_domain(&self.south, 0.0));
        let p10 = self
            .south
            .point_at(remap_unit_to_curve_domain(&self.south, 1.0));
        let p01 = self
            .north
            .point_at(remap_unit_to_curve_domain(&self.north, 0.0));
        let p11 = self
            .north
            .point_at(remap_unit_to_curve_domain(&self.north, 1.0));

        let linear_u = south * (1.0 - v) + north * v;
        let linear_v = west * (1.0 - u) + east * u;
        let bilinear = p00 * ((1.0 - u) * (1.0 - v))
            + p10 * (u * (1.0 - v))
            + p01 * ((1.0 - u) * v)
            + p11 * (u * v);
        linear_u + linear_v - bilinear
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = self.point_at((u + eps).clamp(0.0, 1.0), v)
            - self.point_at((u - eps).clamp(0.0, 1.0), v);
        let dv = self.point_at(u, (v + eps).clamp(0.0, 1.0))
            - self.point_at(u, (v - eps).clamp(0.0, 1.0));
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

/// De Boor's algorithm in homogeneous 4D space.
/// Returns `[wx, wy, wz, w]` (not divided by w yet).
pub(crate) fn de_boor_homo(
    degree: usize,
    knots: &[f64],
    points: &[DVec3],
    weights: &[f64],
    t: f64,
) -> [f64; 4] {
    let n = points.len();
    if n == 0 {
        return [0.0; 4];
    }
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };
    // OCCT BSplCLib: a NULL weight array (non-rational curve) is treated as
    // all weights equal to 1.0. Missing trailing entries are likewise 1.0.
    let mut d: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights.get(idx).copied().unwrap_or(1.0);
            [points[idx].x * w, points[idx].y * w, points[idx].z * w, w]
        })
        .collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }
    d[degree]
}

/// De Boor's algorithm in homogeneous 3D space for 2D rational curves.
/// Returns `[wx, wy, w]` (not divided by w yet).
pub(crate) fn de_boor_homo_2d(
    degree: usize,
    knots: &[f64],
    points: &[DVec2],
    weights: &[f64],
    t: f64,
) -> [f64; 3] {
    let n = points.len();
    if n == 0 {
        return [0.0; 3];
    }
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };
    // OCCT BSplCLib: a NULL weight array (non-rational curve) is treated as
    // all weights equal to 1.0. Missing trailing entries are likewise 1.0.
    let mut d: Vec<[f64; 3]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights.get(idx).copied().unwrap_or(1.0);
            [points[idx].x * w, points[idx].y * w, w]
        })
        .collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }
    d[degree]
}

/// De Boor's algorithm for rational B-spline evaluation.
/// Returns the 3D point at parameter `t`.
pub(crate) fn de_boor(degree: usize, knots: &[f64], points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 {
        return DVec3::ZERO;
    }

    // Find knot span index k such that knots[k] <= t < knots[k+1]
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };

    // Initialize homogeneous control points for the span
    let mut d: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let idx = k - degree + j;
            let idx = idx.min(n - 1);
            let w = weights[idx];
            [points[idx].x * w, points[idx].y * w, points[idx].z * w, w]
        })
        .collect();

    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }

    let w = d[degree][3];
    if w.abs() < 1e-15 {
        DVec3::ZERO
    } else {
        DVec3::new(d[degree][0] / w, d[degree][1] / w, d[degree][2] / w)
    }
}

/// De Boor's algorithm for rational B-spline evaluation in 2D parameter space.
/// Returns the 2D point at parameter `t`. Identical logic to `de_boor` with DVec2.
pub(crate) fn de_boor_2d(degree: usize, knots: &[f64], points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 {
        return DVec2::ZERO;
    }

    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };

    // Homogeneous control points [x*w, y*w, w]
    // OCCT BSplCLib: a NULL weight array (non-rational curve) is treated as
    // all weights equal to 1.0. Missing trailing entries are likewise 1.0.
    let mut d: Vec<[f64; 3]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights.get(idx).copied().unwrap_or(1.0);
            [points[idx].x * w, points[idx].y * w, w]
        })
        .collect();

    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }

    let w = d[degree][2];
    if w.abs() < 1e-15 {
        DVec2::ZERO
    } else {
        DVec2::new(d[degree][0] / w, d[degree][1] / w)
    }
}

/// Analytic tangent for a rational B-Spline curve (NURBS) using the quotient rule.
///
/// The derivative of C(t) = A(t)/W(t) is:
///   C'(t) = (A'(t) - W'(t)*C(t)) / W(t)
///
/// A'(t) and W'(t) are degree-(p-1) B-Splines with control points:
///   A'_i = p * (w_{i+1}*P_{i+1} - w_i*P_i) / (t_{i+p+1} - t_{i+1})
///   W'_i = p * (w_{i+1} - w_i)              / (t_{i+p+1} - t_{i+1})
///
/// Returns the unnormalised derivative vector (caller normalises if needed).
pub(crate) fn bspline_tangent_analytic(
    degree: usize,
    knots: &[f64],
    points: &[DVec3],
    weights: &[f64],
    t: f64,
) -> DVec3 {
    let n = points.len();
    if n < 2 || degree == 0 {
        return DVec3::ZERO;
    }

    // OCCT: a NULL weight array (non-rational curve) is treated as all weights
    // equal to 1.0 (BSplCLib weight accessor, cf. math/bspl.rs wgt).
    let ws: Vec<f64> = if weights.is_empty() { vec![1.0; n] } else { weights.to_vec() };

    let p = degree as f64;
    let m = n - 1; // number of derivative control points

    let mut a_prime: Vec<DVec3> = Vec::with_capacity(m);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(m); // scalar stored in .x
    for i in 0..m {
        let denom = knots[i + degree + 1] - knots[i + 1];
        if denom.abs() < 1e-15 {
            a_prime.push(DVec3::ZERO);
            w_prime.push(DVec3::ZERO);
        } else {
            let s = p / denom;
            a_prime.push(s * (ws[i + 1] * points[i + 1] - ws[i] * points[i]));
            w_prime.push(DVec3::new(s * (ws[i + 1] - ws[i]), 0.0, 0.0));
        }
    }

    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0f64; m];

    // A'(t): non-rational B-Spline of degree p-1
    let a_prime_t = de_boor(degree - 1, deriv_knots, &a_prime, &unit, t);
    // W'(t): scalar B-Spline of degree p-1 (embedded in .x)
    let w_prime_t = de_boor(degree - 1, deriv_knots, &w_prime, &unit, t).x;

    // W(t) and C(t) from the homogeneous evaluation
    let h = crate::math::bspl::de_boor_homo(degree, knots, points, weights, t);
    let w_t = h[3];
    if w_t.abs() < 1e-15 {
        return DVec3::ZERO;
    }
    let c_t = DVec3::new(h[0] / w_t, h[1] / w_t, h[2] / w_t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

pub(crate) fn bspline_tangent_analytic_2d(
    degree: usize,
    knots: &[f64],
    points: &[DVec2],
    weights: &[f64],
    t: f64,
) -> DVec2 {
    let n = points.len();
    if n < 2 || degree == 0 {
        return DVec2::ZERO;
    }

    // OCCT: a NULL weight array (non-rational curve) is treated as all weights
    // equal to 1.0 (BSplCLib weight accessor, cf. math/bspl.rs wgt).
    let ws: Vec<f64> = if weights.is_empty() { vec![1.0; n] } else { weights.to_vec() };

    let p = degree as f64;
    let m = n - 1;

    let mut a_prime = Vec::with_capacity(m);
    let mut w_prime = Vec::with_capacity(m);
    for i in 0..m {
        let denom = knots[i + degree + 1] - knots[i + 1];
        if denom.abs() < 1e-15 {
            a_prime.push(DVec2::ZERO);
            w_prime.push(DVec2::ZERO);
        } else {
            let s = p / denom;
            a_prime.push(s * (ws[i + 1] * points[i + 1] - ws[i] * points[i]));
            w_prime.push(DVec2::new(s * (ws[i + 1] - ws[i]), 0.0));
        }
    }

    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0; m];
    let a_prime_t = de_boor_2d(degree - 1, deriv_knots, &a_prime, &unit, t);
    let w_prime_t = de_boor_2d(degree - 1, deriv_knots, &w_prime, &unit, t).x;

    let h = de_boor_homo_2d(degree, knots, points, weights, t);
    let w_t = h[2];
    if w_t.abs() < 1e-15 {
        return DVec2::ZERO;
    }
    let c_t = DVec2::new(h[0] / w_t, h[1] / w_t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

/// Analytic tangent for a rational Bezier curve using the quotient rule.
///
/// The derivative of a degree-n Bezier is a degree-(n-1) Bezier with:
///   A'_i = n*(w_{i+1}*P_{i+1} - w_i*P_i)
///   W'_i = n*(w_{i+1} - w_i)
pub(crate) fn bezier_tangent_analytic(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n < 2 {
        return DVec3::ZERO;
    }
    let deg = (n - 1) as f64;

    let mut a_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        a_prime.push(deg * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
        w_prime.push(DVec3::new(deg * (weights[i + 1] - weights[i]), 0.0, 0.0));
    }

    let unit = vec![1.0f64; n - 1];
    let a_prime_t = de_casteljau_3d(&a_prime, &unit, t);
    let w_prime_t = de_casteljau_3d(&w_prime, &unit, t).x;

    // W(t): evaluate weights as scalar Bezier (embed in .x with unit weights)
    let w_pts: Vec<DVec3> = weights.iter().map(|&w| DVec3::new(w, 0.0, 0.0)).collect();
    let w_unit = vec![1.0f64; n]; // n elements to match w_pts
    let w_t = de_casteljau_3d(&w_pts, &w_unit, t).x;
    if w_t.abs() < 1e-15 {
        return DVec3::ZERO;
    }

    // C(t) from the standard rational evaluation
    let c_t = de_casteljau_3d(points, weights, t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

/// Analytic derivative for a rational Bezier curve in 2D.
/// Same formula as `bezier_tangent_analytic` but operating on DVec2.
pub(crate) fn bezier_tangent_analytic_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n < 2 {
        return DVec2::ZERO;
    }
    let deg = (n - 1) as f64;
    let mut a_prime: Vec<DVec2> = Vec::with_capacity(n - 1);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        a_prime.push(deg * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
        w_prime.push(DVec3::new(deg * (weights[i + 1] - weights[i]), 0.0, 0.0));
    }
    let unit = vec![1.0f64; n - 1];
    let a_prime_t = de_casteljau_2d(&a_prime, &unit, t);
    // Evaluate w'(t) — use DVec3 to embed w' scalar in .x
    let w_prime_t = de_casteljau_3d(&w_prime, &unit, t).x;
    let w_pts: Vec<DVec3> = weights.iter().map(|&w| DVec3::new(w, 0.0, 0.0)).collect();
    let w_unit = vec![1.0f64; n];
    let w_t = de_casteljau_3d(&w_pts, &w_unit, t).x;
    if w_t.abs() < 1e-15 {
        return DVec2::ZERO;
    }
    let c_t = de_casteljau_2d(points, weights, t);
    a_prime_t - (w_prime_t * c_t) / w_t
}

impl CurveEval for BSplineCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        crate::math::bspl::de_boor(
            self.degree, &self.knots, &self.control_points, &self.weights, t,
        )
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t,
        ).normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t,
        )
    }
    fn derivative2_at(&self, t: f64) -> DVec3 {
        let h = 1e-5;
        (crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t + h,
        ) - crate::math::bspl::bspline_tangent(
            self.degree, &self.knots, &self.control_points, &self.weights, t - h,
        )) / (2.0 * h)
    }
    fn derivative3_at(&self, t: f64) -> DVec3 {
        let h = 1e-4;
        (self.derivative2_at(t + h) - self.derivative2_at(t - h)) / (2.0 * h)
    }
    fn default_domain(&self) -> [f64; 2] {
        let d = self.degree;
        let n = self.knots.len();
        if n < 2 * d + 2 { return [0.0, 1.0]; }
        [self.knots[d], self.knots[n - d - 1]]
    }
}

impl SurfaceEval for BSplineSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        // Tensor product rational evaluation (NURBS):
        // 1. For each v-column, evaluate the u-direction NURBS in homogeneous coords
        //    -> get (wx, wy, wz, w) for each column index.
        // 2. Collect column weights and weighted positions.
        // 3. Run de Boor in v on the homogeneous results, then divide by weight.
        let n_u = self.control_points.len();
        if n_u == 0 {
            return DVec3::ZERO;
        }
        let n_v = self.control_points[0].len();
        if n_v == 0 {
            return DVec3::ZERO;
        }
        // Step 1: evaluate each v-column in the u direction -> homogeneous 4-vector
        let col_homo: Vec<[f64; 4]> = (0..n_v)
            .map(|j| {
                let pts: Vec<DVec3> = (0..n_u).map(|i| self.control_points[i][j]).collect();
                let wts: Vec<f64> = (0..n_u).map(|i| self.weights[i][j]).collect();
                crate::math::bspl::de_boor_homo(self.degree_u, &self.knots_u, &pts, &wts, u)
            })
            .collect();
        // Step 2: build the v-direction "control points" and "weights" from col_homo
        let v_pts: Vec<DVec3> = col_homo
            .iter()
            .map(|h| {
                let w = h[3];
                if w.abs() < 1e-15 {
                    DVec3::ZERO
                } else {
                    DVec3::new(h[0] / w, h[1] / w, h[2] / w)
                }
            })
            .collect();
        let v_wts: Vec<f64> = col_homo.iter().map(|h| h[3]).collect();
        // Step 3: rational de Boor in v
        crate::math::bspl::de_boor(self.degree_v, &self.knots_v, &v_pts, &v_wts, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let [_u0, u1, _v0, v1] = self.default_domain();
        let du = if u + eps <= u1 {
            self.point_at(u + eps, v) - self.point_at(u, v)
        } else {
            self.point_at(u, v) - self.point_at(u - eps, v)
        };
        let dv = if v + eps <= v1 {
            self.point_at(u, v + eps) - self.point_at(u, v)
        } else {
            self.point_at(u, v) - self.point_at(u, v - eps)
        };
        let n = du.cross(dv);
        let len = n.length();
        if len < 1e-15 { DVec3::Z } else { n / len }
    }
    fn default_domain(&self) -> [f64; 4] {
        let du = self.degree_u;
        let dv = self.degree_v;
        let nu = self.knots_u.len();
        let nv = self.knots_v.len();
        let u0 = if nu > du { self.knots_u[du] } else { 0.0 };
        let u1 = if nu > du + 1 {
            self.knots_u[nu - du - 1]
        } else {
            1.0
        };
        let v0 = if nv > dv { self.knots_v[dv] } else { 0.0 };
        let v1 = if nv > dv + 1 {
            self.knots_v[nv - dv - 1]
        } else {
            1.0
        };
        [u0, u1, v0, v1]
    }
}

// --- Curve2dEval implementations ---

pub(crate) fn de_casteljau_3d(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 {
        return DVec3::ZERO;
    }
    // Work in homogeneous coordinates [x*w, y*w, z*w, w]
    let mut d: Vec<[f64; 4]> = points
        .iter()
        .zip(weights)
        .map(|(p, &w)| [p.x * w, p.y * w, p.z * w, w])
        .collect();
    for r in 1..n {
        for j in 0..n - r {
            let next = d[j + 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(next.iter()) {
                *elem = (1.0 - t) * *elem + t * p;
            }
        }
    }
    let w = d[0][3];
    if w.abs() < 1e-15 {
        DVec3::ZERO
    } else {
        DVec3::new(d[0][0] / w, d[0][1] / w, d[0][2] / w)
    }
}

/// De Casteljau algorithm for rational Bezier curve evaluation in 2D.
pub(crate) fn de_casteljau_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 {
        return DVec2::ZERO;
    }
    let mut d: Vec<[f64; 3]> = points
        .iter()
        .zip(weights)
        .map(|(p, &w)| [p.x * w, p.y * w, w])
        .collect();
    for r in 1..n {
        for j in 0..n - r {
            let next = d[j + 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(next.iter()) {
                *elem = (1.0 - t) * *elem + t * p;
            }
        }
    }
    let w = d[0][2];
    if w.abs() < 1e-15 {
        DVec2::ZERO
    } else {
        DVec2::new(d[0][0] / w, d[0][1] / w)
    }
}

impl CurveEval for BezierCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        de_casteljau_3d(&self.control_points, &self.weights, t)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        bezier_tangent_analytic(&self.control_points, &self.weights, t).normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        bezier_tangent_analytic(&self.control_points, &self.weights, t)
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }
}

impl SurfaceEval for BezierSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let n_u = self.control_points.len();
        if n_u == 0 {
            return DVec3::ZERO;
        }
        let n_v = self.control_points[0].len();
        if n_v == 0 {
            return DVec3::ZERO;
        }
        // Apply de Casteljau in u for each v-column, producing n_v intermediate points
        let row_points: Vec<DVec3> = (0..n_v)
            .map(|j| {
                let col_pts: Vec<DVec3> = (0..n_u).map(|i| self.control_points[i][j]).collect();
                let col_wts: Vec<f64> = (0..n_u).map(|i| self.weights[i][j]).collect();
                de_casteljau_3d(&col_pts, &col_wts, u)
            })
            .collect();
        let unit_wts = vec![1.0; n_v];
        de_casteljau_3d(&row_points, &unit_wts, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        let n = du.cross(dv);
        let len = n.length();
        if len < 1e-15 { DVec3::Z } else { n / len }
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, v| acc * v as f64)
}

fn trinomial_coeff(n: usize, i: usize, j: usize, k: usize) -> f64 {
    factorial(n) / (factorial(i) * factorial(j) * factorial(k))
}

impl SurfaceEval for TriBezierSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let degree = self.control_points.len().saturating_sub(1);
        if self.control_points.is_empty() || self.weights.len() != self.control_points.len() {
            return DVec3::ZERO;
        }

        let w = 1.0 - u - v;
        let mut homo = [0.0; 4];
        for (i, row) in self.control_points.iter().enumerate() {
            if row.len() != degree + 1 - i
                || self.weights.get(i).map(|r| r.len()) != Some(row.len())
            {
                return DVec3::ZERO;
            }
            for (j, point) in row.iter().enumerate() {
                let k = degree - i - j;
                let basis = trinomial_coeff(degree, i, j, k)
                    * u.powi(i as i32)
                    * v.powi(j as i32)
                    * w.powi(k as i32);
                let weight = self.weights[i][j];
                homo[0] += basis * weight * point.x;
                homo[1] += basis * weight * point.y;
                homo[2] += basis * weight * point.z;
                homo[3] += basis * weight;
            }
        }

        if homo[3].abs() < 1e-15 {
            DVec3::ZERO
        } else {
            DVec3::new(homo[0] / homo[3], homo[1] / homo[3], homo[2] / homo[3])
        }
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

impl SurfaceEval for RuledSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let start = self.start.point_at(u);
        let end = self.end.point_at(u);
        start.lerp(end, v)
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = self.end.point_at(u) - self.start.point_at(u);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        let [u0, u1] = self.start.default_domain();
        [u0, u1, 0.0, 1.0]
    }
}

impl Curve2dEval for BezierCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        de_casteljau_2d(&self.control_points, &self.weights, t)
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        bezier_tangent_analytic_2d(&self.control_points, &self.weights, t)
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }
}

impl CurveEval for OffsetCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let base_pt = self.basis.point_at(t);
        let tangent = self.basis.tangent_at(t);
        let perp = tangent.cross(self.offset_dir);
        let perp_len = perp.length();
        if perp_len < 1e-15 {
            return base_pt;
        }
        base_pt + self.offset_distance * (perp / perp_len)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let eps = 1e-6;
        let [t0, t1] = self.basis.default_domain();
        let t_lo = (t - eps).max(t0);
        let t_hi = (t + eps).min(t1);
        let dp = self.point_at(t_hi) - self.point_at(t_lo);
        let len = dp.length();
        if len < 1e-15 { DVec3::X } else { dp / len }
    }
    fn default_domain(&self) -> [f64; 2] {
        self.basis.default_domain()
    }
}

impl SurfaceEval for OffsetSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let base_pt = self.basis.point_at(u, v);
        let n = self.basis.normal_at(u, v);
        base_pt + self.offset_distance * n
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        // Offset preserves the normal direction (first-order approximation)
        self.basis.normal_at(u, v)
    }
    fn default_domain(&self) -> [f64; 4] {
        self.basis.default_domain()
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
    fn reversed_parameter(&self, t: f64) -> f64 {
        -t
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
        let minor = DVec2::new(-self.major_dir.y, self.major_dir.x);
        self.center
            + self.major_dir * (self.major_radius * t.cos())
            + minor * (self.minor_radius * t.sin())
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        let minor = DVec2::new(-self.major_dir.y, self.major_dir.x);
        (-self.major_radius * t.sin() * self.major_dir + self.minor_radius * t.cos() * minor)
            .normalize()
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        let minor = DVec2::new(-self.major_dir.y, self.major_dir.x);
        -self.major_radius * t.sin() * self.major_dir + self.minor_radius * t.cos() * minor
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

impl Curve2dEval for BSplineCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        crate::math::bspl::de_boor_2d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        crate::math::bspl::bspline_tangent_2d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        let d = self.degree;
        let n = self.knots.len();
        if n > 2 * d {
            [self.knots[d], self.knots[n - d - 1]]
        } else if n >= 2 {
            [self.knots[0], self.knots[n - 1]]
        } else {
            [0.0, 1.0]
        }
    }
}

impl Curve2dEval for Curve2d {
    fn point_at(&self, t: f64) -> DVec2 {
        match self {
            Curve2d::Trimmed(tc) => tc.point_at(t),
            Curve2d::Line(c) => c.point_at(t),
            Curve2d::Circle(c) => c.point_at(t),
            Curve2d::Ellipse(c) => c.point_at(t),
            Curve2d::CircleInvolute(c) => c.point_at(t),
            Curve2d::Parabola(c) => c.point_at(t),
            Curve2d::Hyperbola(c) => c.point_at(t),
            Curve2d::ArchimedeanSpiral(c) => c.point_at(t),
            Curve2d::LogarithmicSpiral(c) => c.point_at(t),
            Curve2d::SineWave(c) => c.point_at(t),
            Curve2d::BSpline(c) => c.point_at(t),
            Curve2d::Bezier(c) => c.point_at(t),
            Curve2d::Offset(c) => c.point_at(t),
            Curve2d::AHTBezier(c) => c.point_at(t),
            Curve2d::TBezier(c) => c.point_at(t),
        }
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        match self {
            Curve2d::Trimmed(tc) => tc.tangent_at(t),
            Curve2d::Line(c) => c.tangent_at(t),
            Curve2d::Circle(c) => c.tangent_at(t),
            Curve2d::Ellipse(c) => c.tangent_at(t),
            Curve2d::CircleInvolute(c) => c.tangent_at(t),
            Curve2d::Parabola(c) => c.tangent_at(t),
            Curve2d::Hyperbola(c) => c.tangent_at(t),
            Curve2d::ArchimedeanSpiral(c) => c.tangent_at(t),
            Curve2d::LogarithmicSpiral(c) => c.tangent_at(t),
            Curve2d::SineWave(c) => c.tangent_at(t),
            Curve2d::BSpline(c) => c.tangent_at(t),
            Curve2d::Bezier(c) => c.tangent_at(t),
            Curve2d::Offset(c) => c.tangent_at(t),
            Curve2d::AHTBezier(c) => c.derivative_at(t).normalize_or_zero(),
            Curve2d::TBezier(c) => c.derivative_at(t).normalize_or_zero(),
        }
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        match self {
            Curve2d::Trimmed(tc) => tc.derivative_at(t),
            Curve2d::Line(c) => c.derivative_at(t),
            Curve2d::Circle(c) => c.derivative_at(t),
            Curve2d::Ellipse(c) => c.derivative_at(t),
            Curve2d::CircleInvolute(c) => c.derivative_at(t),
            Curve2d::Parabola(c) => c.derivative_at(t),
            Curve2d::Hyperbola(c) => c.derivative_at(t),
            Curve2d::ArchimedeanSpiral(c) => c.derivative_at(t),
            Curve2d::LogarithmicSpiral(c) => c.derivative_at(t),
            Curve2d::SineWave(c) => c.derivative_at(t),
            Curve2d::BSpline(c) => c.derivative_at(t),
            Curve2d::Bezier(c) => c.derivative_at(t),
            Curve2d::Offset(c) => c.derivative_at(t),
            Curve2d::AHTBezier(c) => c.derivative_at(t),
            Curve2d::TBezier(c) => c.derivative_at(t),
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        match self {
            Curve2d::Trimmed(tc) => [tc.t_min, tc.t_max],
            Curve2d::Line(_) => [f64::NEG_INFINITY, f64::INFINITY],
            Curve2d::Circle(_) => [0.0, 2.0 * PI],
            Curve2d::Ellipse(_) => [0.0, 2.0 * PI],
            Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => [f64::NEG_INFINITY, f64::INFINITY],
            Curve2d::CircleInvolute(_) => [0.0, 10.0],
            Curve2d::ArchimedeanSpiral(_) => [0.0, 6.0 * PI],
            Curve2d::LogarithmicSpiral(_) => [0.0, 4.0 * PI],
            Curve2d::SineWave(_) => [-10.0, 10.0],
            Curve2d::BSpline(c) => {
                let d = c.degree;
                let n = c.knots.len();
                if n > 2 * d {
                    [c.knots[d], c.knots[n - d - 1]]
                } else {
                    [0.0, 1.0]
                }
            }
            Curve2d::Bezier(_) => [0.0, 1.0],
            Curve2d::Offset(c) => c.basis.default_domain(),
            Curve2d::AHTBezier(_) => [0.0, 1.0],
            Curve2d::TBezier(c) => [0.0, std::f64::consts::PI / c.alpha],
        }
    }
}

impl Curve2dEval for TrimmedCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        let t_clamped = t.clamp(self.t_min, self.t_max);
        match self.curve.as_ref() {
            Curve2d::BSpline(_) | Curve2d::Bezier(_) => {
                let span = self.t_max - self.t_min;
                if span > 0.0 {
                    let t_norm = (t_clamped - self.t_min) / span;
                    self.curve.point_at(t_norm)
                } else {
                    self.curve.point_at(0.0)
                }
            }
            _ => self.curve.point_at(t_clamped),
        }
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        let t_clamped = t.clamp(self.t_min, self.t_max);
        match self.curve.as_ref() {
            Curve2d::BSpline(_) | Curve2d::Bezier(_) => {
                let span = self.t_max - self.t_min;
                if span > 0.0 {
                    let t_norm = (t_clamped - self.t_min) / span;
                    self.curve.tangent_at(t_norm)
                } else {
                    self.curve.tangent_at(0.0)
                }
            }
            _ => self.curve.tangent_at(t_clamped),
        }
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        let t_clamped = t.clamp(self.t_min, self.t_max);
        self.curve.derivative_at(t_clamped)
    }
    fn default_domain(&self) -> [f64; 2] {
        [self.t_min, self.t_max]
    }
}

impl Curve2dEval for OffsetCurve2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let base_pt = self.basis.point_at(t);
        // Compute tangent via finite differences
        let eps = 1e-6;
        let t_hi = t + eps;
        let t_lo = t - eps;
        let dp = self.basis.point_at(t_hi) - self.basis.point_at(t_lo);
        let tangent = dp.normalize_or_zero();
        // OCCT-aligned right-hand normal: Z_cross_tangent = (Ty, -Tx)
        let normal = DVec2::new(tangent.y, -tangent.x);
        base_pt + self.offset_distance * normal
    }
}

fn aht_basis_values(t: f64, alg_deg: usize, alpha: f64, beta: f64) -> Vec<f64> {
    // Basis: {1, t, ..., t^k, sinh(αt), cosh(αt), sin(βt), cos(βt)}
    let mut basis = Vec::new();
    // Polynomial part: 1, t, t^2, ..., t^k
    let mut tp = 1.0;
    for _ in 0..=alg_deg {
        basis.push(tp);
        tp *= t;
    }
    // Hyperbolic part: sinh(αt), cosh(αt)
    if alpha > 0.0 {
        let a = alpha * t;
        basis.push(a.sinh());
        basis.push(a.cosh());
    }
    // Trigonometric part: sin(βt), cos(βt)
    if beta > 0.0 {
        let b = beta * t;
        basis.push(b.sin());
        basis.push(b.cos());
    }
    basis
}

impl Curve2dEval for AHTBezierCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        let basis = aht_basis_values(t, self.alg_degree, self.alpha, self.beta);
        let n = self.control_points.len().min(basis.len());
        if self.weights.is_empty() {
            // Non-rational: straight sum
            let mut pt = DVec2::ZERO;
            for i in 0..n {
                pt += self.control_points[i] * basis[i];
            }
            pt
        } else {
            // Rational: weighted sum / weight sum
            let mut pt = DVec2::ZERO;
            let mut wsum = 0.0;
            for i in 0..n {
                let w = if i < self.weights.len() {
                    self.weights[i]
                } else {
                    1.0
                };
                pt += self.control_points[i] * (w * basis[i]);
                wsum += w * basis[i];
            }
            if wsum.abs() > 1e-15 { pt / wsum } else { pt }
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }
}

impl Curve2dEval for TBezierCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        // Basis: {1, cos(αt), sin(αt), cos(2αt), sin(2αt), ..., cos(n·αt), sin(n·αt)}
        let n = self.order;
        let at = self.alpha * t;
        let mut pt = DVec2::ZERO;
        let mut wsum = 0.0;
        let has_weights = !self.weights.is_empty();
        // Constant basis = 1
        let w0 = if has_weights { self.weights[0] } else { 1.0 };
        pt += self.control_points[0] * w0;
        wsum += w0;
        for i in 1..=n {
            let fi = i as f64;
            let c = (fi * at).cos();
            let s = (fi * at).sin();
            let idx_c = 2 * i - 1;
            let idx_s = 2 * i;
            if idx_c < self.control_points.len() {
                let wc = if has_weights && idx_c < self.weights.len() {
                    self.weights[idx_c]
                } else {
                    1.0
                };
                pt += self.control_points[idx_c] * (wc * c);
                wsum += wc * c;
            }
            if idx_s < self.control_points.len() {
                let ws = if has_weights && idx_s < self.weights.len() {
                    self.weights[idx_s]
                } else {
                    1.0
                };
                pt += self.control_points[idx_s] * (ws * s);
                wsum += ws * s;
            }
        }
        if has_weights && wsum.abs() > 1e-15 {
            pt / wsum
        } else {
            pt
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, std::f64::consts::PI / self.alpha]
    }
}

// --- Curve2d helper methods ---

impl Curve2d {
    /// Unwrap through a [`Curve2d::Trimmed`] layer, returning a reference to
    /// the innermost curve. If not trimmed, returns `self` unchanged.
    pub fn inner(&self) -> &Curve2d {
        match self {
            Curve2d::Trimmed(tc) => tc.curve.as_ref(),
            other => other,
        }
    }
}

