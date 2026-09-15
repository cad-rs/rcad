//! Enum-level accessors and dispatch impls: the OCCT `IsKind` / `DownCast`
//! equivalents (`is_conic` / `as_conic` / `as_bounded` / ...) on `Curve3`,
//! `Curve2d` and `Surface3`, the trimmed-curve evaluation impls and the
//! per-variant dispatch of `CurveEval for Curve3`, `Curve2dEval for Curve2d`
//! and `SurfaceEval for Surface3`.
//!
//! Extracted verbatim from `geom/eval.rs` (project Rule 5 file-size split);
//! the evaluation traits stay re-exported from `geom/eval.rs`.
use crate::geom::*;
use std::f64::consts::PI;

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
    /// OCCT `Geom_Curve::ReversedParameter` (Geom_Curve.hxx L87, overridden per
    /// concrete curve): the parameter of the same point on the reversed curve.
    /// `Geom_Line::ReversedParameter(U) = -U` (Geom_Line.cxx L163),
    /// `Geom_Circle`/`Geom_Ellipse` = `2*PI - U` (Geom_Circle.cxx L184,
    /// Geom_Ellipse.cxx L199), `Geom_TrimmedCurve` delegates to its basis
    /// (Geom_TrimmedCurve.cxx L88-91), `Geom_BSplineCurve` = `UFirst + ULast - U`
    /// (with the reversed parameter range for a periodic curve).
    fn reversed_parameter(&self, t: f64) -> f64 {
        match self {
            Curve3::Line(c) => c.reversed_parameter(t),
            Curve3::Circle(c) => c.reversed_parameter(t),
            Curve3::Ellipse(c) => c.reversed_parameter(t),
            Curve3::BSpline(c) => c.reversed_parameter(t),
            Curve3::Bezier(c) => c.reversed_parameter(t),
            Curve3::Offset(c) => c.reversed_parameter(t),
            Curve3::Hyperbola(c) => c.reversed_parameter(t),
            Curve3::Parabola(c) => c.reversed_parameter(t),
            Curve3::CircularHelix(c) => c.reversed_parameter(t),
            Curve3::SineWave(c) => c.reversed_parameter(t),
            Curve3::Trimmed(tc) => tc.reversed_parameter(t),
        }
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
            // OCCT Geom_RectangularTrimmedSurface::IsUClosed
            // (Geom_RectangularTrimmedSurface.cxx L559-573): an untrimmed U
            // returns the basis IsUClosed; a trimmed U is closed when the
            // basis is U-periodic and the trim length is a whole number of
            // periods. rcad approximates "untrimmed" by the trim range
            // covering the basis domain.
            Surface3::Trimmed(t) => {
                let basis = t.basis.as_ref();
                let [u1, u2, _, _] = t.trim;
                let [d_u1, d_u2, _, _] = basis.default_domain();
                let full = (u1 - d_u1).abs() < crate::core::precision::PCONFUSION
                    && (u2 - d_u2).abs() < crate::core::precision::PCONFUSION;
                if full {
                    basis.is_u_closed()
                } else if basis.is_u_periodic() {
                    let period = d_u2 - d_u1;
                    let len = u2 - u1;
                    len > crate::core::precision::PCONFUSION
                        && (len % period).abs() <= crate::core::precision::PCONFUSION
                } else {
                    false
                }
            }
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
            // OCCT Geom_RectangularTrimmedSurface::IsVClosed
            // (Geom_RectangularTrimmedSurface.cxx L580-593): same rule as U,
            // with the basis IsVClosed / V-periodicity.
            Surface3::Trimmed(t) => {
                let basis = t.basis.as_ref();
                let [_, _, v1, v2] = t.trim;
                let [_, _, d_v1, d_v2] = basis.default_domain();
                let full = (v1 - d_v1).abs() < crate::core::precision::PCONFUSION
                    && (v2 - d_v2).abs() < crate::core::precision::PCONFUSION;
                if full {
                    basis.is_v_closed()
                } else if basis.is_v_periodic() {
                    let period = d_v2 - d_v1;
                    let len = v2 - v1;
                    len > crate::core::precision::PCONFUSION
                        && (len % period).abs() <= crate::core::precision::PCONFUSION
                } else {
                    false
                }
            }
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
    fn derivative2_at(&self, t: f64) -> DVec2 {
        match self {
            Curve2d::Trimmed(tc) => tc.derivative2_at(t),
            Curve2d::Line(c) => c.derivative2_at(t),
            Curve2d::Circle(c) => c.derivative2_at(t),
            Curve2d::Ellipse(c) => c.derivative2_at(t),
            Curve2d::CircleInvolute(c) => c.derivative2_at(t),
            Curve2d::Parabola(c) => c.derivative2_at(t),
            Curve2d::Hyperbola(c) => c.derivative2_at(t),
            Curve2d::ArchimedeanSpiral(c) => c.derivative2_at(t),
            Curve2d::LogarithmicSpiral(c) => c.derivative2_at(t),
            Curve2d::SineWave(c) => c.derivative2_at(t),
            Curve2d::BSpline(c) => c.derivative2_at(t),
            Curve2d::Bezier(c) => c.derivative2_at(t),
            Curve2d::Offset(c) => c.derivative2_at(t),
            Curve2d::AHTBezier(c) => c.derivative2_at(t),
            Curve2d::TBezier(c) => c.derivative2_at(t),
        }
    }
    fn derivative3_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_Curve::D3 is virtual: forward per variant (the same
        // union dispatch as the 3D Curve3 enum).  Without this arm the
        // offset EvalD2 read the basis D3 through the finite-difference
        // trait default even for the exact BSpline basis.
        match self {
            Curve2d::Trimmed(tc) => tc.derivative3_at(t),
            Curve2d::Line(c) => c.derivative3_at(t),
            Curve2d::Circle(c) => c.derivative3_at(t),
            Curve2d::Ellipse(c) => c.derivative3_at(t),
            Curve2d::CircleInvolute(c) => c.derivative3_at(t),
            Curve2d::Parabola(c) => c.derivative3_at(t),
            Curve2d::Hyperbola(c) => c.derivative3_at(t),
            Curve2d::ArchimedeanSpiral(c) => c.derivative3_at(t),
            Curve2d::LogarithmicSpiral(c) => c.derivative3_at(t),
            Curve2d::SineWave(c) => c.derivative3_at(t),
            Curve2d::BSpline(c) => c.derivative3_at(t),
            Curve2d::Bezier(c) => c.derivative3_at(t),
            Curve2d::Offset(c) => c.derivative3_at(t),
            Curve2d::AHTBezier(c) => c.derivative3_at(t),
            Curve2d::TBezier(c) => c.derivative3_at(t),
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        match self {
            Curve2d::Trimmed(tc) => [tc.t_min, tc.t_max],
            // OCCT virtual dispatch: Geom2d_Line / Geom2d_Parabola /
            // Geom2d_Hyperbola FirstParameter/LastParameter return
            // -/+Precision::Infinite() (Geom2d_Line.cxx L144/L151,
            // Geom2d_Parabola.cxx L130/L137, Geom2d_Hyperbola.cxx L147/L154).
            Curve2d::Line(c) => c.default_domain(),
            Curve2d::Circle(_) => [0.0, 2.0 * PI],
            Curve2d::Ellipse(_) => [0.0, 2.0 * PI],
            Curve2d::Parabola(c) => c.default_domain(),
            Curve2d::Hyperbola(c) => c.default_domain(),
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
    fn reversed_parameter(&self, t: f64) -> f64 {
        // OCCT Geom2d_Curve::ReversedParameter:
        // Line/Parabola/Hyperbola/Bezier/BSpline: -U
        // Circle/Ellipse (periodic): Period - U
        match self {
            Curve2d::Line(c) => c.reversed_parameter(t),
            Curve2d::Circle(c) => c.reversed_parameter(t),
            Curve2d::Ellipse(c) => c.reversed_parameter(t),
            Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => -t,
            // OCCT Geom2d_BSplineCurve::ReversedParameter = -U (non-periodic)
            Curve2d::BSpline(_) => -t,
            Curve2d::Bezier(_) => -t,
            Curve2d::Offset(c) => c.basis.reversed_parameter(t),
            _ => -t,
        }
    }
    fn is_closed(&self) -> bool {
        // OCCT Geom2d_Curve::IsClosed delegated per variant (the enum impl
        // previously fell back to the trait default `false` for every
        // variant, including the periodic Circle/Ellipse).
        match self {
            Curve2d::Trimmed(tc) => tc.is_closed(),
            Curve2d::Line(c) => c.is_closed(),
            Curve2d::Circle(c) => c.is_closed(),
            Curve2d::Ellipse(c) => c.is_closed(),
            Curve2d::Parabola(c) => c.is_closed(),
            Curve2d::Hyperbola(c) => c.is_closed(),
            Curve2d::CircleInvolute(c) => c.is_closed(),
            Curve2d::ArchimedeanSpiral(c) => c.is_closed(),
            Curve2d::LogarithmicSpiral(c) => c.is_closed(),
            Curve2d::SineWave(c) => c.is_closed(),
            Curve2d::BSpline(c) => c.is_closed(),
            Curve2d::Bezier(c) => c.is_closed(),
            Curve2d::Offset(c) => c.basis.is_closed(),
            Curve2d::AHTBezier(c) => c.is_closed(),
            Curve2d::TBezier(c) => c.is_closed(),
        }
    }
    fn is_periodic(&self) -> bool {
        // OCCT Geom2d_Curve::IsPeriodic delegated per variant.
        match self {
            Curve2d::Trimmed(tc) => tc.is_periodic(),
            Curve2d::Line(c) => c.is_periodic(),
            Curve2d::Circle(c) => c.is_periodic(),
            Curve2d::Ellipse(c) => c.is_periodic(),
            Curve2d::Parabola(c) => c.is_periodic(),
            Curve2d::Hyperbola(c) => c.is_periodic(),
            Curve2d::CircleInvolute(c) => c.is_periodic(),
            Curve2d::ArchimedeanSpiral(c) => c.is_periodic(),
            Curve2d::LogarithmicSpiral(c) => c.is_periodic(),
            Curve2d::SineWave(c) => c.is_periodic(),
            Curve2d::BSpline(c) => c.is_periodic(),
            Curve2d::Bezier(c) => c.is_periodic(),
            Curve2d::Offset(c) => c.basis.is_periodic(),
            Curve2d::AHTBezier(c) => c.is_periodic(),
            Curve2d::TBezier(c) => c.is_periodic(),
        }
    }
}

impl Curve2dEval for TrimmedCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        let t_clamped = t.clamp(self.t_min, self.t_max);
        // OCCT Geom2d_TrimmedCurve::Value (Geom2d_TrimmedCurve.cxx): clamp the
        // parameter to [FirstParameter, LastParameter] and delegate to the
        // basis curve — no re-normalization.
        self.curve.point_at(t_clamped)
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        let t_clamped = t.clamp(self.t_min, self.t_max);
        self.curve.tangent_at(t_clamped)
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        let t_clamped = t.clamp(self.t_min, self.t_max);
        self.curve.derivative_at(t_clamped)
    }
    fn default_domain(&self) -> [f64; 2] {
        [self.t_min, self.t_max]
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
