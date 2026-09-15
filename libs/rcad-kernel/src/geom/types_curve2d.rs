//! 2D (parameter-space) curve payload types: line, conics, involute, spirals,
//! sine wave, BSpline / Bezier, trimmed, offset, AHT Bezier and T-Bezier, plus
//! the [`Curve2d`] enum.
//!
//! Extracted verbatim from `geom/mod.rs` (project Rule 5 file-size split);
//! all public paths stay stable through the `pub use` re-exports in `mod.rs`.
use glam::DVec2;
use serde::{Deserialize, Serialize};

use super::{Point2, Vec2};
use crate::geom::eval::bspline_tangent_analytic_2d;

/// A line in 2D parameter space: point + direction (OCCT-aligned: gp_Lin2d / Geom2d_Line).
///
/// OCCT: gp_Dir2d is ALWAYS a unit vector — direction must be normalized.
/// Use `Line2d::new()` which enforces unit direction.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Line2d {
    pub origin: Point2,
    /// Unit direction vector (OCCT gp_Dir2d invariant — always unit length).
    pub direction: Vec2,
}

impl Line2d {
    /// OCCT-aligned: construct from point + direction (direction is normalized to unit).
    pub fn new(origin: DVec2, direction: DVec2) -> Self {
        Line2d {
            origin,
            direction: direction.normalize_or_zero(),
        }
    }

    /// OCCT-aligned: Geom2d_Line::Distance(gp_Pnt2d)
    /// Perpendicular distance from a point to the infinite line.
    pub fn distance(&self, point: DVec2) -> f64 {
        let d = point - self.origin;
        // In 2D: cross = |(p - o) × dir| where |dir| = 1 (unit invariant)
        (d.x * self.direction.y - d.y * self.direction.x).abs()
    }

    /// OCCT-aligned: Geom2d_Line::ReversedParameter(t) → -t
    pub fn reversed_parameter(&self, t: f64) -> f64 {
        -t
    }

    /// OCCT-aligned: Geom2d_Line::IsClosed() → false
    pub fn is_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: Geom2d_Line::IsPeriodic() → false
    pub fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT-aligned: Geom2d_Line::SetDirection(gp_Dir2d)
    pub fn with_direction(&self, direction: DVec2) -> Self {
        Line2d::new(self.origin, direction)
    }

    /// OCCT-aligned: Geom2d_Line::SetLocation(gp_Pnt2d)
    pub fn with_origin(&self, origin: DVec2) -> Self {
        Line2d {
            origin,
            direction: self.direction,
        }
    }

    /// OCCT-aligned: Geom2d_Line::Transform(gp_Trsf2d) — translation
    pub fn translate(&self, offset: DVec2) -> Self {
        Line2d {
            origin: self.origin + offset,
            direction: self.direction,
        }
    }

    /// OCCT-aligned: Geom2d_Line::Transform(gp_Trsf2d) — rotation around center
    pub fn rotate(&self, center: DVec2, angle_rad: f64) -> Self {
        let cos_a = angle_rad.cos();
        let sin_a = angle_rad.sin();
        let p = self.origin - center;
        let origin = center + DVec2::new(p.x * cos_a - p.y * sin_a, p.x * sin_a + p.y * cos_a);
        let dir = DVec2::new(
            self.direction.x * cos_a - self.direction.y * sin_a,
            self.direction.x * sin_a + self.direction.y * cos_a,
        );
        Line2d::new(origin, dir)
    }
}

/// A circle in 2D parameter space.
///
/// OCCT-aligned: gp_Circ2d / Geom2d_Circle.
/// Parametric form: `P(t) = center + x_dir * R*cos(t) + y_dir * R*sin(t)`
/// where `x_dir` and `y_dir` are orthogonal unit vectors defining the
/// orientation frame.  Rotating this frame around `center` by angle `dU`
/// is equivalent to shifting the parameter `t → t + dU`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Circle2d {
    pub center: Point2,
    /// X-axis direction of the local frame (unit vector, cos(0) direction).
    /// Default: (1, 0).
    #[serde(default = "x_dir_default")]
    pub x_dir: Vec2,
    /// Y-axis direction of the local frame (unit vector, sin(0) direction).
    /// Default: (0, 1).
    #[serde(default = "y_dir_default")]
    pub y_dir: Vec2,
    pub radius: f64,
}

fn x_dir_default() -> Vec2 {
    DVec2::X
}
fn y_dir_default() -> Vec2 {
    DVec2::Y
}

impl Circle2d {
    /// Create a circle with identity frame (X=(1,0), Y=(0,1)).
    pub fn new(center: Point2, radius: f64) -> Self {
        Self {
            center,
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius,
        }
    }

    /// Rotate the circle's local frame around its center by `angle` radians.
    /// Equivalent to OCCT `gp_Trsf2d::SetRotation(Center, angle)` followed by
    /// `Geom2d_Circle::Transform(Trsf)`.  After rotation,
    /// `new_curve(t) = old_curve(t + angle)`.
    pub fn rotate_center(&mut self, angle: f64) {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let x = self.x_dir;
        let y = self.y_dir;
        self.x_dir = DVec2::new(x.x * cos_a + y.x * sin_a, x.y * cos_a + y.y * sin_a);
        self.y_dir = DVec2::new(-x.x * sin_a + y.x * cos_a, -x.y * sin_a + y.y * cos_a);
    }
}

/// An ellipse in 2D parameter space.
///
/// Analogous to OCCT `Geom2d_Ellipse`. Used as a PCurve when an edge traces
/// an elliptical path on the parameter domain of an adjacent surface.
///
/// Parametric form: `center + major_dir * a*cos(t) + minor_dir * b*sin(t)`
/// where `minor_dir` is the stored Y direction of the positioning 2D axis
/// (OCCT gp_Ax22d keeps both directions; `gp_Elips2d::Reverse` negates the
/// Y direction and keeps X).  Default domain: `[0, 2π]`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Ellipse2d {
    pub center: Point2,
    pub major_dir: Vec2,
    pub minor_dir: Vec2,
    pub major_radius: f64,
    pub minor_radius: f64,
}

/// A 2D parabola in parameter space.
///
/// OCCT-aligned: gp_Parab2d. Parameterization in local frame:
///   X(t) = t²/(2*p), Y(t) = t
/// where p = focal_param (distance from focus to directrix).
/// Default domain: (-inf, +inf).
/// The parabola is positioned by `origin` (apex), `axis_dir` (symmetry axis),
/// and `focal_param > 0`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Parabola2d {
    pub origin: Point2,
    /// Normalized axis direction (from apex toward focus).
    pub axis_dir: Vec2,
    /// Focal parameter p (> 0). The focus is at distance p/2 from apex along axis_dir.
    pub focal_param: f64,
}

/// A 2D hyperbola branch in parameter space.
///
/// OCCT-aligned: gp_Hypr2d. The branch on the positive side of the major axis.
/// Implicit: X²/a² - Y²/b² = 1 in local frame.
/// Parametric:  X(t) = a*cosh(t), Y(t) = b*sinh(t)
/// Default domain: (-inf, +inf).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Hyperbola2d {
    pub center: Point2,
    /// Normalized major-axis direction.
    pub major_dir: Vec2,
    pub semi_major: f64, // a (transverse semi-axis, > 0)
    pub semi_minor: f64, // b (conjugate semi-axis, > 0)
}

/// A 2D involute of a base circle in parameter space.
///
/// Parametric form around the local x-axis:
/// `x(t) = r * (cos t + t sin t)`
/// `y(t) = r * (sin t - t cos t)`
///
/// The local frame is then rotated by `start_angle` and translated by `center`.
/// This curve is commonly used for gear-tooth flank profiles.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CircleInvolute2d {
    pub center: Point2,
    pub base_radius: f64,
    /// Rotation of the local involute frame in radians.
    pub start_angle: f64,
}

/// A 2D Archimedean spiral in parameter space.
///
/// `r(t) = a + b*t`, `theta(t) = start_angle + t`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ArchimedeanSpiral2d {
    pub center: Point2,
    pub a: f64,
    pub b: f64,
    pub start_angle: f64,
}

/// A 2D logarithmic spiral in parameter space.
///
/// `r(t) = a * exp(b*t)`, `theta(t) = start_angle + t`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct LogarithmicSpiral2d {
    pub center: Point2,
    pub a: f64,
    pub b: f64,
    pub start_angle: f64,
}

/// A 2D sine-wave curve in parameter space.
///
/// Parametric form:
/// `x(t) = t`
/// `y(t) = amplitude * sin(frequency * t + phase)`
///
/// Useful for procedural sketching and for matching OCCT's sine-wave evaluator
/// family in a lightweight form.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SineWave2d {
    pub amplitude: f64,
    pub frequency: f64,
    pub phase: f64,
}

/// A non-uniform rational B-spline curve in 2D parameter space.
///
/// Analogous to OCCT `Geom2d_BSplineCurve`. Used for PCurves: the image of
/// a 3D edge in the (u, v) domain of an adjacent surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSplineCurve2 {
    pub degree: usize,
    /// Full knot vector (with multiplicities expanded).
    pub knots: Vec<f64>,
    pub control_points: Vec<DVec2>,
    /// Homogeneous weights; 1.0 for non-rational.
    pub weights: Vec<f64>,
    /// OCCT Geom2d_BSplineCurve::IsPeriodic() — true for a periodic (unclamped)
    /// B-spline whose first/last knot multiplicity equals the degree and whose
    /// poles wrap with the period. When true, the effective parameter range is
    /// [knots[0], knots[len-1]] (FirstParameter = Knots(1), LastParameter =
    /// Knots(NbKnots)) and evaluations wrap across the seam.
    /// OCCT: Geom2d_BSplineCurve.hxx myPeriodic / IsPeriodic().
    #[serde(default)]
    pub is_periodic: bool,
}

impl BSplineCurve2 {
    /// Returns the unnormalized first derivative at parameter `t`.
    pub fn derivative_at(&self, t: f64) -> DVec2 {
        if self.is_periodic {
            // OCCT Geom2d_BSplineCurve::D1 -> EvalD1 (Geom2d_BSplineCurve_1.cxx
            // L199-226): the periodic PeriodicNormalization/LocateParameter
            // wrap — the analytic stencil has no pole wrap.
            return crate::geom::bspline2d_dn::eval_d1(self, t).d1;
        }
        bspline_tangent_analytic_2d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }

    /// Approximate a sequence of ordered 2D points with a cubic BSpline.
    /// Uses chord-length parameterization and natural end conditions.
    /// The resulting curve has degree 3 and passes through all input points.
    pub fn approximate(points: &[DVec2]) -> Self {
        let n = points.len();
        if n < 2 {
            return Self::degenerate();
        }
        if n == 2 {
            // Linear BSpline
            let knots = vec![0.0, 0.0, 1.0, 1.0];
            return BSplineCurve2 {
                degree: 1,
                knots,
                control_points: vec![points[0], points[1]],
                weights: vec![1.0, 1.0],
                is_periodic: false,
            };
        }
        if n == 3 {
            // Quadratic BSpline
            let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
            return BSplineCurve2 {
                degree: 2,
                knots,
                control_points: points.to_vec(),
                weights: vec![1.0; n],
                is_periodic: false,
            };
        }

        // Cubic BSpline interpolation with chord-length parameterization
        let degree = 3;

        // Compute chord-length parameters
        let mut params = vec![0.0_f64; n];
        for i in 1..n {
            let d = (points[i] - points[i - 1]).length();
            params[i] = params[i - 1] + d.max(1e-15);
        }
        let total = params[n - 1];
        for p in &mut params {
            *p /= total;
        }

        // Clamped knot vector with multiplicity = degree+1 at ends
        let n_knots = n + degree + 1;
        let mut knots = vec![0.0_f64; n_knots];

        // First degree+1 knots = 0
        for k in &mut knots[..=degree] {
            *k = params[0];
        }

        // Interior knots (averaging of params)
        for j in 1..n - degree {
            let mut sum = 0.0;
            for i in j..j + degree {
                sum += params[i];
            }
            knots[j + degree] = sum / (degree as f64);
        }

        // Last degree+1 knots = 1
        for k in &mut knots[n_knots - degree - 1..] {
            *k = params[n - 1];
        }

        BSplineCurve2 {
            degree,
            knots,
            control_points: points.to_vec(),
            weights: vec![1.0; n],
            is_periodic: false,
        }
    }

    /// Approximate a closed sequence of 2D points with a periodic cubic BSpline.
    pub fn approximate_closed(points: &[DVec2]) -> Self {
        if points.is_empty() {
            return Self::degenerate();
        }
        // Add the start point at the end to close the loop
        let mut closed_pts = points.to_vec();
        if points.len() > 2 && (points[0] - points[points.len() - 1]).length() > 1e-15 {
            closed_pts.push(points[0]);
        }
        Self::approximate(&closed_pts)
    }

    fn degenerate() -> Self {
        BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::ZERO, DVec2::ZERO],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        }
    }
}

/// A rational or non-rational Bezier curve in 2D parameter space.
///
/// Analogous to OCCT `Geom2d_BezierCurve`. Domain is `[0, 1]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BezierCurve2 {
    pub control_points: Vec<DVec2>,
    /// Homogeneous weights; 1.0 for non-rational.
    pub weights: Vec<f64>,
}

/// A trimmed 2D curve — wraps a curve with a restricted parameter range.
///
/// Analogous to OCCT `Geom2d_TrimmedCurve`. `point_at(t)` clamps `t` to
/// `[t_min, t_max]` before delegating, emulating OCCT's behavior of returning
/// the endpoint value for out-of-range parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimmedCurve2 {
    /// The underlying curve (full-range expression).
    pub curve: Box<Curve2d>,
    /// Lower bound of the valid parameter range.
    pub t_min: f64,
    /// Upper bound of the valid parameter range.
    pub t_max: f64,
}

/// A 2D curve offset from a base curve by a fixed distance along the right-hand normal.
///
/// `P(t) = P_base(t) + offset_distance * N(t)`
/// where `N(t) = Z_cross_T(t) = (Ty, -Tx)` is the unit normal pointing to the
/// right of the direction of travel (OCCT convention). The tangent `T(t)` is
/// computed via finite differences when the base curve does not provide an
/// analytic derivative.
///
/// Analogous to OCCT `Geom2d_OffsetCurve`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffsetCurve2d {
    /// The basis curve.
    pub basis: Box<Curve2d>,
    /// Offset distance (positive = right of travel direction, OCCT convention).
    pub offset_distance: f64,
}

/// A 2D Algebraic-Hyperbolic-Trigonometric (AHT) Bezier curve.
///
/// Basis functions: `{1, t, ..., t^k, sinh(α·t), cosh(α·t), sin(β·t), cos(β·t)}`
/// where `k = alg_degree`. OCCT equivalent: `Geom2dEval_AHTBezierCurve`.
///
/// Number of poles = `alg_degree + 1 + (alpha > 0 ? 2 : 0) + (beta > 0 ? 2 : 0)`.
/// Domain is `[0, 1]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AHTBezierCurve2 {
    pub control_points: Vec<DVec2>,
    /// Homogeneous weights; empty for non-rational.
    pub weights: Vec<f64>,
    /// Algebraic degree k (polynomial part `{1, t, ..., t^k}`).
    pub alg_degree: usize,
    /// Hyperbolic coefficient α (0 = no hyperbolic terms).
    pub alpha: f64,
    /// Trigonometric coefficient β (0 = no trigonometric terms).
    pub beta: f64,
}

/// A 2D Trigonometric Bezier (T-Bezier) curve.
///
/// Basis functions: `{1, cos(t), sin(t), cos(2·t), sin(2·t), ..., cos(n·t), sin(n·t)}`
/// where `n = order`. OCCT equivalent: `Geom2dEval_TBezierCurve`.
///
/// Number of poles = `2·order + 1`. Domain is `[0, π/α]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TBezierCurve2 {
    pub control_points: Vec<DVec2>,
    /// Homogeneous weights; empty for non-rational.
    pub weights: Vec<f64>,
    /// Order n (trigonometric degree).
    pub order: usize,
    /// Frequency-scaling factor α (> 0). Domain = `[0, π/α]`.
    pub alpha: f64,
}

/// A curve defined in the 2D parameter space (u, v) of a surface.
///
/// Used for PCurves: the image of a 3D edge on the parameter domain of an
/// adjacent face surface. Analogous to OCCT `Geom2d_Curve`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Curve2d {
    Line(Line2d),
    Circle(Circle2d),
    Ellipse(Ellipse2d),
    Parabola(Parabola2d),
    Hyperbola(Hyperbola2d),
    CircleInvolute(CircleInvolute2d),
    ArchimedeanSpiral(ArchimedeanSpiral2d),
    LogarithmicSpiral(LogarithmicSpiral2d),
    SineWave(SineWave2d),
    BSpline(BSplineCurve2),
    Bezier(BezierCurve2),
    /// Trimmed curve: restricts evaluation to `[t_min, t_max]`.
    /// See [`TrimmedCurve2`] for details.
    Trimmed(TrimmedCurve2),
    /// A 2D curve offset from a base curve by a fixed distance along the left normal.
    ///
    /// `P(t) = P_base(t) + offset_distance * N(t)`
    /// where `N(t) = Rot90(T(t)) = (-Ty, Tx)` is the unit normal pointing to the
    /// left of the direction of travel.
    ///
    /// Analogous to OCCT `Geom2d_OffsetCurve`.
    Offset(OffsetCurve2d),
    /// Algebraic-Hyperbolic-Trigonometric Bezier curve (AHT Bezier).
    /// OCCT equivalent: `Geom2dEval_AHTBezierCurve`.
    AHTBezier(AHTBezierCurve2),
    /// Trigonometric Bezier curve (T-Bezier).
    /// OCCT equivalent: `Geom2dEval_TBezierCurve`.
    TBezier(TBezierCurve2),
}
