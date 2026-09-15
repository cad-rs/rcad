//! 3D curve payload types: line, conics, helix, sine wave, BSpline / Bezier,
//! trimmed and offset curves, plus the [`Curve3`] enum.
//!
//! Extracted verbatim from `geom/mod.rs` (project Rule 5 file-size split);
//! all public paths stay stable through the `pub use` re-exports in `mod.rs`.
use glam::DVec3;
use serde::{Deserialize, Serialize};

use super::{Point3, Vec3};
use crate::geom::eval::bspline_tangent_analytic;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Line3 {
    pub origin: Point3,
    /// Unit direction vector (OCCT gp_Dir invariant — always unit length).
    /// Use Line3::new() for normalized construction.
    pub direction: Vec3,
}

impl Line3 {
    /// OCCT-aligned: construct from point and direction (direction normalized to unit).
    pub fn new(origin: DVec3, direction: DVec3) -> Self {
        Line3 {
            origin,
            direction: direction.normalize_or_zero(),
        }
    }

    /// OCCT-aligned: gp_Lin::Distance(gp_Pnt) — perpendicular distance.
    pub fn distance(&self, point: DVec3) -> f64 {
        let d = point - self.origin;
        // |d × direction| / |direction| — direction is unit, so denominator = 1
        d.cross(self.direction).length()
    }

    /// OCCT-aligned: Geom_Line::ReversedParameter(t) = -t
    pub fn reversed_parameter(&self, t: f64) -> f64 {
        -t
    }

    /// OCCT-aligned: Geom_Line::IsClosed() = false
    pub fn is_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: Geom_Line::IsPeriodic() = false
    pub fn is_periodic(&self) -> bool {
        false
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Circle3 {
    pub center: Point3,
    pub normal: Vec3,
    #[serde(default = "circle3_x_dir_default")]
    pub x_dir: Vec3,
    #[serde(default = "circle3_y_dir_default")]
    pub y_dir: Vec3,
    pub radius: f64,
}

fn stable_x_dir(normal: Vec3) -> Vec3 {
    // OCCT-aligned: gp_Ax2 reference direction for circle's local frame.
    // Use X as reference; if N is parallel to X (|N·X| ≥ 1-1e-12), use Z instead.
    // u_dir = (ref - N * (N·ref)).normalize()
    let ref_dir = if normal.x.abs() > 1.0 - 1e-12 {
        DVec3::Z
    } else {
        DVec3::X
    };
    (ref_dir - normal * ref_dir.dot(normal)).normalize_or_zero()
}
fn circle3_x_dir_default() -> Vec3 {
    DVec3::X
}
fn circle3_y_dir_default() -> Vec3 {
    DVec3::Y
}

impl Circle3 {
    /// OCCT-aligned: construct a circle with orthonormal frame.
    pub fn new(center: Point3, normal: Vec3, radius: f64) -> Self {
        let normal = normal.normalize_or_zero();
        let x_dir = stable_x_dir(normal);
        let y_dir = normal.cross(x_dir).normalize();
        Self {
            center,
            normal,
            x_dir,
            y_dir,
            radius,
        }
    }

    /// OCCT-aligned: gp_Circ::Distance(gp_Pnt) — min distance from point to circle curve.
    pub fn distance(&self, point: DVec3) -> f64 {
        let d = point - self.center;
        let axis_dist = d.dot(self.normal);
        let planar = d - axis_dist * self.normal;
        let planar_dist = planar.length();
        let radial_diff = (planar_dist - self.radius).abs();
        (axis_dist * axis_dist + radial_diff * radial_diff).sqrt()
    }

    /// OCCT gp_Circ(gp_Ax2(P, N, X)) — construct a circle with an explicit
    /// reference direction for the local frame (x_dir = ref_dir projected onto
    /// the circle plane). Mirrors BRepPrim_OneAxis::TopEdge/BottomEdge which
    /// build the cap circles with Axes().XDirection() as the reference.
    pub fn new_with_ref_dir(center: Point3, normal: Vec3, radius: f64, ref_dir: Vec3) -> Self {
        let normal = normal.normalize_or_zero();
        let ref_rej = ref_dir - normal * ref_dir.dot(normal);
        let x_dir = if ref_rej.length_squared() < 1e-12 {
            stable_x_dir(normal)
        } else {
            ref_rej.normalize()
        };
        let y_dir = normal.cross(x_dir).normalize();
        Self {
            center,
            normal,
            x_dir,
            y_dir,
            radius,
        }
    }

    pub fn rotate_frame(&mut self, angle: f64) {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let x = self.x_dir;
        let y = self.y_dir;
        self.x_dir = DVec3::new(
            x.x * cos_a + y.x * sin_a,
            x.y * cos_a + y.y * sin_a,
            x.z * cos_a + y.z * sin_a,
        );
        self.y_dir = DVec3::new(
            -x.x * sin_a + y.x * cos_a,
            -x.y * sin_a + y.y * cos_a,
            -x.z * sin_a + y.z * cos_a,
        );
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Ellipse3 {
    pub center: Point3,
    pub normal: Vec3,
    pub major_dir: Vec3,
    pub major_radius: f64,
    pub minor_radius: f64,
}

/// A non-uniform rational B-spline curve in 3D.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSplineCurve3 {
    pub degree: usize,
    /// Full knot vector (with multiplicities expanded).
    pub knots: Vec<f64>,
    pub control_points: Vec<DVec3>,
    /// Homogeneous weights; 1.0 for non-rational.
    pub weights: Vec<f64>,
    /// OCCT Geom_BSplineCurve::IsPeriodic() — true for a periodic (unclamped)
    /// B-spline whose first/last knot multiplicity equals the degree and whose
    /// poles wrap with the period. When true, the effective parameter range is
    /// [knots[degree], knots[n-degree-1]] and evaluations wrap across the seam.
    /// OCCT: Geom_BSplineCurve.hxx myPeriodic / IsPeriodic().
    #[serde(default)]
    pub is_periodic: bool,
}

impl BSplineCurve3 {
    /// Returns the unnormalized first derivative at parameter `t`.
    pub fn derivative_at(&self, t: f64) -> DVec3 {
        bspline_tangent_analytic(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }

    /// ✅ OCCT-aligned: C2-continuous knot intervals.
    /// Returns the knot values bounding each C2-continuous span.
    /// Equivalent to OCCT's NbIntervals(curve, GeomAbs_C2).
    /// Between consecutive returned values the curve has C2 continuity.
    pub fn c2_intervals(&self) -> Vec<f64> {
        let d = self.degree;
        let n = self.knots.len();
        let t_min = self.knots[d];
        let t_max = self.knots[n - d - 1];
        let mut boundaries = Vec::new();
        boundaries.push(t_min);
        // Skip the first knot multiplicity (the first d knots are clamped)
        let mut i = d + 1;
        while i < n - d {
            let k = self.knots[i];
            // Count multiplicity at this knot
            let mut m = 1_usize;
            while i + 1 < n - d && (self.knots[i + 1] - k).abs() < 1e-15 {
                i += 1;
                m += 1;
            }
            // OCCT: C2 boundary when multiplicity < degree
            if m < d && k > t_min && k < t_max {
                boundaries.push(k);
            }
            i += 1;
        }
        boundaries.push(t_max);
        // Deduplicate near-equal boundaries
        boundaries.dedup_by(|a, b| (*a - *b).abs() < 1e-14);
        boundaries
    }
}

/// A rational or non-rational Bezier curve in 3D.
///
/// Evaluated via de Casteljau's algorithm. Domain is always `[0.0, 1.0]`.
/// Analogous to OCCT `Geom_BezierCurve`.
///
/// Note: a Bezier curve of degree n is equivalent to a B-spline of degree n
/// with knot vector `[0, ..., 0, 1, ..., 1]` (n+1 times each).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BezierCurve3 {
    pub control_points: Vec<DVec3>,
    /// Homogeneous weights; 1.0 for non-rational (polynomial Bezier).
    pub weights: Vec<f64>,
}

/// A 3D hyperbola defined by center, normal, semi-transverse axis `a`, and
/// semi-conjugate axis `b`.  Parametric form:
///
///   P(t) = center + a * cosh(t) * major_dir + b * sinh(t) * minor_dir
///
/// where `minor_dir = normal × major_dir`.  Domain is `(-∞, +∞)`;
/// the principal branch (t ≥ 0) is on the `+major_dir` side.
/// Analogous to OCCT `Geom_Hyperbola`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Hyperbola3 {
    pub center: Point3,
    pub normal: Vec3,
    pub major_dir: Vec3,
    pub semi_major: f64, // a  (transverse semi-axis)
    pub semi_minor: f64, // b  (conjugate semi-axis)
}

/// A 3D parabola defined by its vertex, axis, and focal parameter `p`
/// (where the focus is at distance `p/2` from the vertex along the axis).
///
///   P(t) = vertex + (t²/(2p)) * axis_dir + t * dir_perp
///
/// where `dir_perp = normal × axis_dir` is the cross-axis direction.
/// Domain is `(-∞, +∞)`.  Analogous to OCCT `Geom_Parabola`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Parabola3 {
    pub vertex: Point3,
    pub normal: Vec3,
    pub axis_dir: Vec3,   // direction from vertex toward focus
    pub focal_param: f64, // p  (= 2 × focal_length)
}

/// A circular helix curve around an axis.
///
/// Parameterization:
/// `P(t) = origin + radius*(cos t * x_axis + sin t * y_axis) + (pitch/(2*pi))*t * axis`
///
/// Analogous to OCCT TKHelix circular helix primitives.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CircularHelix3 {
    pub origin: Point3,
    pub axis: Vec3,
    /// A reference direction orthogonalized against `axis` at evaluation time.
    pub ref_dir: Vec3,
    pub radius: f64,
    /// Axial advance per full revolution (2*pi in parameter).
    pub pitch: f64,
}

/// A 3D sine-wave curve traveling along a baseline direction with amplitude
/// in a perpendicular `amplitude_dir`.
///
/// Parameterization:
/// `P(t) = origin + t * baseline_dir + amplitude * sin(frequency * t + phase) * amplitude_dir`
///
/// `baseline_dir` and `amplitude_dir` should be orthogonal unit vectors.
/// Analogous to OCCT `GeomEval_SineWaveCurve`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SineWave3 {
    pub origin: Point3,
    /// Unit direction along which the parameter `t` advances.
    pub baseline_dir: Vec3,
    /// Unit direction of the sine-wave displacement (orthogonal to `baseline_dir`).
    pub amplitude_dir: Vec3,
    pub amplitude: f64,
    pub frequency: f64,
    pub phase: f64,
}

/// OCCT Geom_TrimmedCurve equivalent: wraps a base Curve3 with parameter domain [first, last].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimmedCurve3 {
    pub curve: Box<Curve3>,
    pub first: f64,
    pub last: f64,
}

impl TrimmedCurve3 {
    pub fn new(curve: Curve3, first: f64, last: f64) -> Self {
        Self {
            curve: Box::new(curve),
            first,
            last,
        }
    }
    pub fn map_param(&self, t: f64) -> f64 {
        t
    }
    pub fn basis_curve(&self) -> &Curve3 {
        &self.curve
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Curve3 {
    Line(Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
    BSpline(BSplineCurve3),
    Bezier(BezierCurve3),
    Offset(OffsetCurve3),
    Hyperbola(Hyperbola3),
    Parabola(Parabola3),
    CircularHelix(CircularHelix3),
    SineWave(SineWave3),
    /// OCCT Geom_TrimmedCurve: a curve bounded to a parameter range [first, last].
    Trimmed(TrimmedCurve3),
}

/// A curve offset from a base curve by a fixed distance in a reference plane.
///
/// `S(t) = basis.point_at(t) + offset_distance * (tangent(t) × offset_dir).normalize()`
///
/// Analogous to OCCT `Geom_OffsetCurve`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffsetCurve3 {
    pub basis: Box<Curve3>,
    /// Offset distance (positive = outward from the curve's "left" side).
    pub offset_distance: f64,
    /// Fixed reference direction (normal to the offset plane).
    /// The offset direction at each point is `(tangent × offset_dir).normalize()`.
    pub offset_dir: Vec3,
}
