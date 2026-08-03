use glam::{DVec2, DVec3};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// Math helpers defined in eval.rs, called by inherent impls in this file
use crate::geom::eval::bspline_tangent_analytic;
use crate::geom::eval::bspline_tangent_analytic_2d;

pub type Point3 = DVec3;
pub type Vec3 = DVec3;
pub type Point2 = DVec2;
pub type Vec2 = DVec2;

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

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Plane {
    pub origin: Point3,
    /// Unit normal vector (OCCT gp_Dir invariant). Use Plane::new() for normalized construction.
    pub normal: Vec3,
    /// U-axis direction (X direction of OCCT gp_Ax3). Orthogonal to normal.
    /// Determines the U=0 direction in the plane's (u, v) parameterization.
    pub u_dir: Vec3,
    /// V-axis direction (Y direction of OCCT gp_Ax3, = normal × u_dir).
    /// Orthogonal to both normal and u_dir.
    pub v_dir: Vec3,
}

impl Plane {
    /// OCCT-aligned: construct from origin and normal.
    /// Equivalent to `gp_Pln(gp_Pnt, gp_Dir)` which internally creates
    /// a `gp_Ax3(P, V)` whose X direction is the unit vector perpendicular
    /// to V having a zero in the coordinate of the smallest |component| of V
    /// (gp_Ax3.cxx L29-80):
    ///   1. A,B,C = V.X,V.Y,V.Z; Aabs,Babs,Cabs = |A|,|B|,|C|
    ///   2. If |B| is smallest:  D = |A|>|C| ? (-C,0,A) : (C,0,-A)
    ///   3. elif |A| is smallest: D = |B|>|C| ? (0,-C,B) : (0,C,-B)
    ///   4. else:                 D = |A|>|B| ? (-B,A,0) : (B,-A,0)
    ///   5. u_dir = D.normalize(); v_dir = N × u_dir
    pub fn new(origin: DVec3, normal: DVec3) -> Self {
        let normal = normal.normalize_or_zero();
        let (a, b, c) = (normal.x, normal.y, normal.z);
        let (aabs, babs, cabs) = (a.abs(), b.abs(), c.abs());
        let d = if babs <= aabs && babs <= cabs {
            if aabs > cabs {
                DVec3::new(-c, 0.0, a)
            } else {
                DVec3::new(c, 0.0, -a)
            }
        } else if aabs <= babs && aabs <= cabs {
            if babs > cabs {
                DVec3::new(0.0, -c, b)
            } else {
                DVec3::new(0.0, c, -b)
            }
        } else if aabs > babs {
            DVec3::new(-b, a, 0.0)
        } else {
            DVec3::new(b, -a, 0.0)
        };
        let u_dir = d.normalize_or_zero();
        let v_dir = normal.cross(u_dir).normalize_or_zero();
        Plane {
            origin,
            normal,
            u_dir,
            v_dir,
        }
    }

    /// OCCT-aligned: construct from origin, normal, and explicit u_dir.
    /// Equivalent to `gp_Pln(gp_Ax3(origin, normal, u_dir))`.
    /// v_dir is computed as `normal × u_dir` (right-handed orthonormal frame).
    pub fn with_axes(origin: DVec3, normal: DVec3, u_dir: DVec3) -> Self {
        let normal = normal.normalize_or_zero();
        let u_dir = u_dir.normalize();
        let v_dir = normal.cross(u_dir).normalize_or_zero();
        Plane {
            origin,
            normal,
            u_dir,
            v_dir,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CylindricalSurface {
    pub origin: Point3,
    pub axis: Vec3,
    pub radius: f64,
    /// Reference direction for u=0 (perpendicular to axis).
    /// Preserved through rotation so UV mapping stays consistent.
    pub ref_dir: Vec3,
}

impl CylindricalSurface {
    /// Create a cylinder with [`any_perpendicular(axis)`](any_perpendicular) as the reference direction.
    pub fn new(origin: Point3, axis: Vec3, radius: f64) -> Self {
        Self {
            origin,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            ref_dir: any_perpendicular(axis),
        }
    }

    /// Create a cylinder with an explicit reference direction for u=0.
    pub fn new_with_ref_dir(origin: Point3, axis: Vec3, radius: f64, ref_dir: Vec3) -> Self {
        Self {
            origin,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            ref_dir: ref_dir.normalize_or_zero(),
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SphericalSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    /// Reference direction for u=0 (perpendicular to axis).
    /// Preserved through rotation so UV mapping stays consistent.
    pub ref_dir: Vec3,
}

impl SphericalSurface {
    /// Direction perpendicular to ref_dir in the equatorial plane (axis × ref_dir).
    /// Used for UV mapping: U = atan2(dot(ref_dir_perp), dot(ref_dir)).
    pub fn ref_dir_perp(&self) -> Vec3 {
        self.axis.cross(self.ref_dir).normalize_or_zero()
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ConicalSurface {
    /// Point on the cone axis where the surface radius equals `radius`.
    ///
    /// Historically this field was used as an apex for zero-radius primitive
    /// cones. For general conical surfaces, the true apex is derived from this
    /// reference point, `radius`, and `half_angle_rad`.
    pub apex: Point3,
    pub axis: Vec3,
    /// Radius of the reference circle at `apex`.
    pub radius: f64,
    pub half_angle_rad: f64,
    /// Reference direction for u=0 (perpendicular to axis).
    ///
    /// OCCT `gp_Cone` carries a full `gp_Ax3` (Position + XDirection); the
    /// XDirection defines the u=0 generatrix.  Preserved through rotation so UV
    /// mapping stays consistent — same convention as [`CylindricalSurface`] and
    /// [`SphericalSurface`].
    #[serde(default = "default_cone_ref_dir")]
    pub ref_dir: Vec3,
}

/// Serde fallback for the cone `ref_dir` field when reading legacy data that
/// predates the field.  New data always carries the explicit reference.
fn default_cone_ref_dir() -> DVec3 {
    DVec3::X
}

impl ConicalSurface {
    /// Create a cone with the default reference direction for u=0
    /// (`any_perpendicular(axis)`, matching OCCT `gp_Ax2` default).
    pub fn new(apex: Point3, axis: Vec3, radius: f64, half_angle_rad: f64) -> Self {
        Self {
            apex,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            half_angle_rad,
            ref_dir: any_perpendicular(axis),
        }
    }

    /// Create a cone with an explicit reference direction for u=0.
    pub fn new_with_ref_dir(apex: Point3, axis: Vec3, radius: f64, half_angle_rad: f64, ref_dir: Vec3) -> Self {
        Self {
            apex,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            half_angle_rad,
            ref_dir: ref_dir.normalize_or_zero(),
        }
    }

    pub fn axis_dir(&self) -> DVec3 {
        self.axis.normalize_or_zero()
    }

    pub fn apex_point(&self) -> DVec3 {
        let tan_half = self.half_angle_rad.tan();
        if tan_half.abs() < 1e-12 {
            self.apex
        } else {
            self.apex - self.axis_dir() * (self.radius / tan_half)
        }
    }

    pub fn axial_from_slant(&self, slant: f64) -> f64 {
        slant * self.half_angle_rad.cos()
    }

    pub fn slant_from_axial(&self, axial: f64) -> f64 {
        let cos_half = self.half_angle_rad.cos();
        if cos_half.abs() < 1e-12 {
            0.0
        } else {
            axial / cos_half
        }
    }

    pub fn radius_at_slant(&self, slant: f64) -> f64 {
        self.radius + slant * self.half_angle_rad.sin()
    }

    pub fn radius_at_axial(&self, axial: f64) -> f64 {
        self.radius + axial * self.half_angle_rad.tan()
    }

    /// UV coordinates of world point `p` relative to this conical surface.
    ///
    /// `u` = azimuth (−π, π], `v` = slant distance from the reference circle at
    /// `self.apex`, matching [`SurfaceEval::point_at`].  When `p` is off the
    /// surface the returned `(u, v)` corresponds to the closest point on the cone.
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let local = p - self.apex;
        let along = local.dot(axis);
        let perp = local - axis * along;
        let radial = perp.length();

        let u = if radial < 1e-15 {
            0.0
        } else {
            let perp_n = perp / radial;
            perp_n.dot(y_ax).atan2(perp_n.dot(x_ax))
        };

        let cos_half = self.half_angle_rad.cos();
        let sin_half = self.half_angle_rad.sin();
        let v = along * cos_half + (radial - self.radius) * sin_half;

        DVec2::new(u, v)
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ToroidalSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub major_radius: f64,
    pub minor_radius: f64,
}

/// An ellipsoidal surface aligned to a local orthonormal frame.
///
/// Parameterization matches sphere-like angles:
/// - `u` = longitude `[0, 2π]`
/// - `v` = colatitude `[0, π]` (0 at +axis pole)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct EllipsoidalSurface {
    pub center: Point3,
    pub axis: Vec3,
    /// Reference direction used to derive the local X axis.
    pub ref_dir: Vec3,
    pub radius_x: f64,
    pub radius_y: f64,
    pub radius_z: f64,
}

/// A classical helicoid surface around an axis.
///
/// Parameterization:
/// `S(u, v) = origin + v * (cos(u) * x_axis + sin(u) * y_axis) + (pitch/(2*pi))*u * axis`
///
/// `u` is the azimuth / screw parameter and `v` is the signed radial distance.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct HelicoidSurface {
    pub origin: Point3,
    pub axis: Vec3,
    /// Reference direction used to derive the local X axis.
    pub ref_dir: Vec3,
    /// Axial advance per full revolution.
    pub pitch: f64,
}

/// A circular pipe/tube surface around a spine curve.
///
/// `u` is the azimuth angle around the local section frame and `v` follows the
/// natural parameter of the spine curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeSurface {
    pub spine: Box<Curve3>,
    /// Initial/reference direction projected onto the normal plane of the
    /// spine tangent at evaluation time.
    pub ref_dir: Vec3,
    pub radius: f64,
}

/// A non-uniform rational B-spline surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSplineSurface {
    pub degree_u: usize,
    pub degree_v: usize,
    /// Full knot vector for u (with multiplicities expanded).
    pub knots_u: Vec<f64>,
    /// Full knot vector for v (with multiplicities expanded).
    pub knots_v: Vec<f64>,
    /// Control point grid [u_index][v_index].
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight grid [u_index][v_index]; 1.0 for non-rational.
    pub weights: Vec<Vec<f64>>,
}

/// Returns `true` if the BSpline surface is planar (degree ≤ 1 in both directions
/// and all control points lie within `tol` of a single plane).
pub fn bspline_is_planar(bsp: &BSplineSurface, tol: f64) -> bool {
    if bsp.control_points.is_empty() {
        return false;
    }

    // Collect all unique control points
    let pts: Vec<DVec3> = bsp
        .control_points
        .iter()
        .flat_map(|row| row.iter())
        .copied()
        .collect();
    if pts.len() < 3 {
        return true; // trivially planar
    }

    // Find the first 3 non-collinear points to define the plane
    let origin = pts[0];
    let mut normal = DVec3::ZERO;
    for i in 1..pts.len() - 1 {
        let d1 = pts[i] - origin;
        let d2 = pts[i + 1] - origin;
        let n = d1.cross(d2);
        if n.length_squared() > tol * tol {
            normal = n.normalize();
            break;
        }
    }
    if normal.length_squared() < 0.5 {
        // All points are collinear — any plane containing the line works
        return true;
    }

    // Check all points lie within tol of the plane
    pts.iter().all(|&p| {
        let d = (p - origin).dot(normal);
        d.abs() <= tol
    })
}

/// Convert a planar BSpline surface to the best-fit `Plane`.
/// The BSpline must satisfy `bspline_is_planar`. Returns the plane from the
/// first 3 non-collinear control points (or any plane if all points collinear).
pub fn bspline_to_plane(bsp: &BSplineSurface) -> Plane {
    let pts: Vec<DVec3> = bsp
        .control_points
        .iter()
        .flat_map(|row| row.iter())
        .copied()
        .collect();

    let origin = pts[0];
    let mut normal = DVec3::Z;
    for i in 1..pts.len() - 1 {
        let d1 = pts[i] - origin;
        let d2 = pts[i + 1] - origin;
        let n = d1.cross(d2);
        if n.length_squared() > 1e-30 {
            normal = n.normalize();
            break;
        }
    }

    Plane::new(origin, normal)
}

/// A rational or non-rational Bezier surface (tensor-product bicubic patch).
///
/// Evaluated by applying de Casteljau in u, then in v. Domain is `[0, 1] × [0, 1]`.
/// Analogous to OCCT `Geom_BezierSurface`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BezierSurface {
    /// Control point grid [u_count][v_count].
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight grid [u_count][v_count]; 1.0 for non-rational.
    pub weights: Vec<Vec<f64>>,
}

/// A triangular rational Bezier surface using barycentric coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriBezierSurface {
    /// Triangular control net rows. Row `i` has `degree + 1 - i` points.
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight rows with the same triangular layout as `control_points`.
    pub weights: Vec<Vec<f64>>,
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

/// A surface offset from a base surface by a fixed distance along the normal.
///
/// `S(u,v) = basis.point_at(u,v) + offset_distance * basis.normal_at(u,v)`
///
/// The offset normal is the same as the basis normal. Domain equals the basis domain.
/// Analogous to OCCT `Geom_OffsetSurface`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffsetSurface {
    pub basis: Box<Surface3>,
    /// Offset distance along the outward normal (positive = outward).
    pub offset_distance: f64,
}

/// A rectangular trimmed surface: a base surface restricted to the UV box
/// `[u1, u2] × [v1, v2]`.
///
/// Evaluation delegates fully to the basis surface; only the reported domain
/// changes. Analogous to OCCT `Geom_RectangularTrimmedSurface`.
///
/// Appears in STEP as `RECTANGULAR_TRIMMED_SURFACE`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimmedSurface {
    /// The underlying surface being trimmed.
    pub basis: Box<Surface3>,
    /// Trim bounds `[u1, u2, v1, v2]`.
    pub trim: [f64; 4],
}

impl TrimmedSurface {
    pub fn new(basis: Surface3, u1: f64, u2: f64, v1: f64, v2: f64) -> Self {
        Self {
            basis: Box::new(basis),
            trim: [u1, u2, v1, v2],
        }
    }
}

/// Surface formed by translating a 3D profile curve along a direction.
/// S(u,v) = profile.point_at(u) + v * direction
/// Analogous to OCCT Geom_SurfaceOfLinearExtrusion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearExtrusionSurface {
    pub profile: Box<Curve3>,
    /// Normalized extrusion direction.
    pub direction: Vec3,
}

/// Surface formed by rotating a 3D profile curve around an axis.
/// S(u,v) = rotate(profile.point_at(v), axis_origin, axis_dir, angle=u)
/// Analogous to OCCT Geom_SurfaceOfRevolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RevolutionSurface {
    pub profile: Box<Curve3>,
    pub axis_origin: Point3,
    /// Normalized rotation axis direction.
    pub axis_dir: Vec3,
}

/// Surface linearly interpolating between two 3D curves with a shared parameter domain.
/// S(u,v) = lerp(start.point_at(u), end.point_at(u), v)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RuledSurface {
    pub start: Box<Curve3>,
    pub end: Box<Curve3>,
}

/// A Coons patch blending four boundary curves over `[0,1] x [0,1]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoonsSurface {
    /// Boundary curve at `v = 0`, parameterized by `u`.
    pub south: Box<Curve3>,
    /// Boundary curve at `v = 1`, parameterized by `u`.
    pub north: Box<Curve3>,
    /// Boundary curve at `u = 0`, parameterized by `v`.
    pub west: Box<Curve3>,
    /// Boundary curve at `u = 1`, parameterized by `v`.
    pub east: Box<Curve3>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Surface3 {
    Plane(Plane),
    Cylinder(CylindricalSurface),
    Sphere(SphericalSurface),
    Cone(ConicalSurface),
    Torus(ToroidalSurface),
    Ellipsoid(EllipsoidalSurface),
    Helicoid(HelicoidSurface),
    Pipe(PipeSurface),
    BSpline(BSplineSurface),
    LinearExtrusion(LinearExtrusionSurface),
    Revolution(RevolutionSurface),
    Ruled(RuledSurface),
    Coons(CoonsSurface),
    Bezier(BezierSurface),
    TriBezier(TriBezierSurface),
    Offset(OffsetSurface),
    Trimmed(TrimmedSurface),
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum PrimitiveSolid {
    Box {
        width: f64,
        height: f64,
        depth: f64,
    },
    Sphere {
        radius: f64,
    },
    Cylinder {
        radius: f64,
        height: f64,
    },
    Cone {
        base_radius: f64,
        height: f64,
    },
    Torus {
        major_radius: f64,
        minor_radius: f64,
    },
}

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
/// where `minor_dir = rotate_ccw_90(major_dir)`.  Default domain: `[0, 2π]`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Ellipse2d {
    pub center: Point2,
    pub major_dir: Vec2,
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
}

impl BSplineCurve2 {
    /// Returns the unnormalized first derivative at parameter `t`.
    pub fn derivative_at(&self, t: f64) -> DVec2 {
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

/// OCCT-aligned: returns a vector perpendicular to `v` using the gp_Ax2
/// reference direction convention (project X onto plane; fallback to Z).
/// Stable for any non-zero input.
pub fn any_perpendicular(v: DVec3) -> DVec3 {
    // OCCT gp_Ax2 reference direction selection:
    // Use X (1,0,0) as default; if |v·X| ≥ 1-1e-12 (v parallel to X), use Z (0,0,1).
    let ref_dir = if v.x.abs() > 1.0 - 1e-12 {
        DVec3::Z
    } else {
        DVec3::X
    };
    // Project reference onto plane perpendicular to v: ref - v*(v·ref)
    let perp = ref_dir - v * ref_dir.dot(v);
    perp.normalize_or_zero()
}

/// OCCT-aligned: 90° counter-clockwise rotation in 2D (gp_Dir2d::Rotated).
/// Equivalent to cross product with (0, 0, 1) in 2D homogeneous form.
pub fn turn_2d(v: DVec2) -> DVec2 {
    DVec2::new(-v.y, v.x)
}

fn orthonormal_frame(axis: DVec3, ref_dir: DVec3) -> (DVec3, DVec3, DVec3) {
    let axis = axis.normalize_or_zero();
    let mut x_axis = ref_dir - axis * ref_dir.dot(axis);
    if x_axis.length_squared() <= 1e-24 {
        x_axis = any_perpendicular(axis);
    } else {
        x_axis = x_axis.normalize();
    }
    let y_axis = axis.cross(x_axis).normalize_or_zero();
    (axis, x_axis, y_axis)
}

/// Parametric evaluation of a 3D curve: `t -> Point3`.
///
/// Mirrors OCCT `Geom_Curve::Value(t)` / `D1(t)`.
pub trait CurveEval {
    /// Point on the curve at parameter `t`.
    fn point_at(&self, t: f64) -> DVec3;
    /// Unit tangent vector at parameter `t`.
    fn tangent_at(&self, t: f64) -> DVec3;
    /// First derivative (non-unit velocity vector) at parameter `t`.
    /// OCCT-aligned: Extrema_CurveTool::D1 — used by Extrema_LocateExtPC
    /// for the Newton method solving g(u) = (P-C)·C' = 0.
    /// Default: 6-point central difference (accurate to O(h⁴)).
    fn derivative_at(&self, t: f64) -> DVec3 {
        // 6-point stencil: f'(t) ≈ [f(t-2h)-8f(t-h)+8f(t+h)-f(t+2h)] / (12h)
        let h = 1e-6;
        let fp2 = self.point_at(t + 2.0 * h);
        let fp1 = self.point_at(t + h);
        let fm1 = self.point_at(t - h);
        let fm2 = self.point_at(t - 2.0 * h);
        (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)
    }
    /// Natural parameter domain `[t_min, t_max]`.
    /// Lines use `[NEG_INFINITY, INFINITY]`; circles/ellipses use `[0, 2π]`.
    fn default_domain(&self) -> [f64; 2];

    /// OCCT-aligned: IsClosed — true for periodic curves where start == end (circle, ellipse).
    fn is_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: IsPeriodic — true for curves with cyclic parameter (circle, ellipse).
    fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT-aligned: ReversedParameter(t) — parameter of the same geometric point
    /// when traversed in the opposite direction.
    /// For periodic curves (circle): `period - t` (mod period).
    /// For non-periodic curves: `-t` (line).
    fn reversed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: D2(t) — second derivative d²P/dt² at parameter `t`.
    ///
    /// Default: 5-point central difference of `point_at`:
    /// ```text
    /// f''(t) = [-f(t+2h) + 16f(t+h) - 30f(t) + 16f(t-h) - f(t-2h)] / (12h²)
    /// ```
    /// Override with analytic formula when available (Line3 → 0, Circle3 → -R·N, etc.).
    fn derivative2_at(&self, t: f64) -> DVec3 {
        let h = 1e-4;
        let fp2 = self.point_at(t + 2.0 * h);
        let fp1 = self.point_at(t + h);
        let f = self.point_at(t);
        let fm1 = self.point_at(t - h);
        let fm2 = self.point_at(t - 2.0 * h);
        (-fp2 + 16.0 * fp1 - 30.0 * f + 16.0 * fm1 - fm2) / (12.0 * h * h)
    }

    /// OCCT-aligned: D3(t) — third derivative d³P/dt³ at parameter `t`.
    ///
    /// Default: central difference of `derivative2_at`.
    fn derivative3_at(&self, t: f64) -> DVec3 {
        let h = 1e-4;
        (self.derivative2_at(t + h) - self.derivative2_at(t - h)) / (2.0 * h)
    }

    /// Signed curvature at parameter `t`.
    ///
    /// For 3D curves: `k = |r' × r''| / |r'|³`.
    /// Returns 0 when the velocity is zero (degenerate point).
    /// OCCT-aligned: computed via `D1 × D2 / |D1|³`.
    fn curvature_at(&self, t: f64) -> f64 {
        let d1 = self.derivative_at(t);
        let d2 = self.derivative2_at(t);
        let speed = d1.length();
        if speed < 1e-15 {
            return 0.0;
        }
        d1.cross(d2).length() / (speed * speed * speed)
    }

    /// OCCT-aligned: `Geom_Curve::TransformedParameter(t, T)`.
    ///
    /// Returns the parameter on the transformed curve corresponding to parameter
    /// `t` on the original curve after transformation `T`. Used for curve-on-surface
    /// evaluation after a `TopLoc_Location` transform.
    ///
    /// Default: identity (`t` unchanged — correct for isometric transformations).
    /// For scaling transformations, override with `t / scale_factor`.
    fn transformed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: `Geom_Curve::ParametricTransformation(T)`.
    ///
    /// Returns the scale factor for parametric transformation. Used when computing
    /// parametric tolerance after a `TopLoc_Location` transformation.
    ///
    /// Default: 1.0 (correct for isometric / uniform scaling).
    fn parametric_transformation(&self) -> f64 {
        1.0
    }
}

/// Parametric evaluation of a 3D surface: `(u, v) -> Point3`.
///
/// Mirrors OCCT `Geom_Surface::Value(u, v)`.
pub trait SurfaceEval {
    /// Point on the surface at parameter `(u, v)`.
    fn point_at(&self, u: f64, v: f64) -> DVec3;
    /// Outward unit normal at parameter `(u, v)`.
    fn normal_at(&self, u: f64, v: f64) -> DVec3;
    /// Natural parameter domain `[u_min, u_max, v_min, v_max]`.
    fn default_domain(&self) -> [f64; 4];
    /// First partial derivatives `(point, dP/du, dP/dv)` at `(u, v)`.
    /// Default: finite-difference approximation (2-point, 1e-6 step).
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let eps = 1e-6;
        let p = self.point_at(u, v);
        let pu = self.point_at(u + eps, v);
        let pv = self.point_at(u, v + eps);
        (p, (pu - p) / eps, (pv - p) / eps)
    }

    /// OCCT-aligned: IsUClosed / IsVClosed — true if surface is closed in that direction.
    fn is_u_closed(&self) -> bool {
        false
    }
    fn is_v_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: IsUPeriodic / IsVPeriodic — true if parameter is cyclic.
    fn is_u_periodic(&self) -> bool {
        false
    }
    fn is_v_periodic(&self) -> bool {
        false
    }

    /// OCCT-aligned: UReversedParameter / VReversedParameter — parameter in reverse direction.
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        t
    }
    fn v_reversed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: D2(u,v) — second-order partial derivatives.
    ///
    /// Returns `(P, dP/du, dP/dv, d²P/du², d²P/dudv, d²P/dv²)`.
    ///
    /// Default: 3-point central difference for each second-order term:
    /// ```text
    /// Puu  = [P(u+h,v) - 2P(u,v) + P(u-h,v)] / h²
    /// Pvv  = [P(u,v+h) - 2P(u,v) + P(u,v-h)] / h²
    /// Puv  = [P(u+h,v+h) - P(u+h,v-h) - P(u-h,v+h) + P(u-h,v-h)] / (4h²)
    /// ```
    /// Override with analytic formula for known surface types.
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let h = 1e-5;
        let (p, pu, pv) = self.derivatives(u, v);
        let p_up = self.point_at(u + h, v);
        let p_um = self.point_at(u - h, v);
        let p_vp = self.point_at(u, v + h);
        let p_vm = self.point_at(u, v - h);
        let p_pp = self.point_at(u + h, v + h);
        let p_pm = self.point_at(u + h, v - h);
        let p_mp = self.point_at(u - h, v + h);
        let p_mm = self.point_at(u - h, v - h);
        let h2 = h * h;
        let puu = (p_up - 2.0 * p + p_um) / h2;
        let pvv = (p_vp - 2.0 * p + p_vm) / h2;
        let puv = (p_pp - p_pm - p_mp + p_mm) / (4.0 * h2);
        (p, pu, pv, puu, puv, pvv)
    }
}

/// Parametric evaluation of a 2D curve (PCurve): `t -> Point2`.
///
/// ✅ OCCT-aligned: Geom2d_Curve (provides Value/D0/D1/D2/D3 + domain).
pub trait Curve2dEval {
    /// Point on the 2D curve at parameter `t` (OCCT `Value(t)` / `D0(t)`).
    fn point_at(&self, t: f64) -> DVec2;

    /// Unit tangent vector at parameter `t` (OCCT `D1(t).normalize()`).
    /// Default: finite-difference approximation.
    fn tangent_at(&self, t: f64) -> DVec2 {
        let eps = 1e-7;
        let dp = self.point_at(t + eps) - self.point_at(t - eps);
        dp.normalize_or_zero()
    }

    /// First derivative (velocity) vector at parameter `t` (OCCT `D1(t)`).
    /// Default: central-difference approximation.
    fn derivative_at(&self, t: f64) -> DVec2 {
        let eps = 1e-7;
        (self.point_at(t + eps) - self.point_at(t - eps)) / (2.0 * eps)
    }

    /// Natural parameter domain `[t_min, t_max]` (OCCT `FirstParameter() / LastParameter()`).
    fn default_domain(&self) -> [f64; 2] {
        [f64::NEG_INFINITY, f64::INFINITY]
    }

    /// OCCT-aligned: IsClosed — true for closed curves (circle, ellipse).
    fn is_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: IsPeriodic — true for periodic curves (circle, ellipse).
    fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT-aligned: ReversedParameter(t) — parameter of the same point in reverse.
    fn reversed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: D2(t) — second derivative d²P/dt² at parameter `t`.
    /// Default: 5-point central difference of `point_at`.
    fn derivative2_at(&self, t: f64) -> DVec2 {
        let h = 1e-4;
        let fp2 = self.point_at(t + 2.0 * h);
        let fp1 = self.point_at(t + h);
        let f = self.point_at(t);
        let fm1 = self.point_at(t - h);
        let fm2 = self.point_at(t - 2.0 * h);
        (-fp2 + 16.0 * fp1 - 30.0 * f + 16.0 * fm1 - fm2) / (12.0 * h * h)
    }

    /// OCCT-aligned: D3(t) — third derivative at parameter `t`.
    /// Default: central difference of `derivative2_at`.
    fn derivative3_at(&self, t: f64) -> DVec2 {
        let h = 1e-4;
        (self.derivative2_at(t + h) - self.derivative2_at(t - h)) / (2.0 * h)
    }

    /// Signed curvature at parameter `t` in the 2D plane.
    ///
    /// `k = (x'y'' - y'x'') / (x'² + y'²)^(3/2)`.
    /// Positive = counter-clockwise turning. Returns 0 at degenerate points.
    fn curvature_at(&self, t: f64) -> f64 {
        let d1 = self.derivative_at(t);
        let d2 = self.derivative2_at(t);
        let speed_sq = d1.length_squared();
        if speed_sq < 1e-30 {
            return 0.0;
        }
        (d1.x * d2.y - d1.y * d2.x) / (speed_sq * speed_sq.sqrt())
    }
}

/// OCCT-aligned: `Geom_Conic` intermediate abstract class.
///
/// Groups all conic curves (Circle, Ellipse, Hyperbola, Parabola) with their
/// shared geometric properties: position frame, eccentricity, and local axes.
///
/// In OCCT boolean algorithms this grouping is used for dynamic type checks
/// like `IsKind(STANDARD_TYPE(Geom_Conic))` before calling conic-specific
/// methods such as `XAxis()` / `YAxis()`.
pub trait ConicEval: CurveEval {
    /// The reference point of the conic:
    ///   Circle/Ellipse/Hyperbola → center
    ///   Parabola → vertex
    fn position(&self) -> DVec3;

    /// The normal of the conic's plane (gp_Ax2::Direction).
    fn normal(&self) -> DVec3;

    /// OCCT-aligned: eccentricity of the conic.
    ///   Circle → 0
    ///   Ellipse → sqrt(1 - (b/a)²)  (0 < e < 1)
    ///   Hyperbola → sqrt(1 + (b/a)²) (e > 1)
    ///   Parabola → 1
    fn eccentricity(&self) -> f64;

    /// OCCT-aligned: XAxis() — the local X direction (major axis for Circle/Ellipse).
    fn x_axis(&self) -> DVec3;

    /// OCCT-aligned: YAxis() — the local Y direction (minor axis for Circle/Ellipse).
    fn y_axis(&self) -> DVec3;
}

/// OCCT-aligned: `Geom_BoundedCurve` intermediate abstract class.
///
/// Groups bounded curves (BSpline, Bezier) that always have a finite domain.
/// All `BoundedCurve` types have `default_domain()` returning finite values.
pub trait BoundedCurveEval: CurveEval {
    /// The degree of the underlying polynomial/rational representation.
    fn degree(&self) -> usize;
}

// --- ConicEval implementations ---

pub trait Conic2dEval: Curve2dEval {
    /// Center or vertex position.
    fn position(&self) -> DVec2;
    /// OCCT-aligned: eccentricity.
    fn eccentricity(&self) -> f64;
    /// Local X-axis direction (major axis for Circle/Ellipse).
    fn x_axis(&self) -> DVec2;
    /// Local Y-axis direction (minor axis for Circle/Ellipse).
    fn y_axis(&self) -> DVec2;
}

/// OCCT-aligned: `Geom2d_BoundedCurve` intermediate abstract class.
pub trait BoundedCurve2dEval: Curve2dEval {
    fn degree(&self) -> usize;
}

pub trait ElementarySurfaceEval: SurfaceEval {
    /// Origin of the local coordinate system (gp_Ax3::Location).
    fn position(&self) -> DVec3;
    /// Normal / axis direction (gp_Ax3::Direction).
    fn axis_dir(&self) -> DVec3;
    /// U-direction (gp_Ax3::XDirection).
    fn x_axis(&self) -> DVec3;
    /// V-direction (gp_Ax3::YDirection).
    fn y_axis(&self) -> DVec3;
}

/// OCCT-aligned: `Geom_BoundedSurface` intermediate abstract class.
///
/// Groups bounded surfaces (BSpline, Bezier) whose domain is always finite.
pub trait BoundedSurfaceEval: SurfaceEval {
    fn degree_u(&self) -> usize;
    fn degree_v(&self) -> usize;
}

/// OCCT-aligned: `Geom_SweptSurface` intermediate abstract class.
///
/// Groups swept surfaces (LinearExtrusion, Revolution) that share a profile curve.
pub trait SweptSurfaceEval: SurfaceEval {
    fn profile(&self) -> &Curve3;
}

pub fn transform_curve(curve: &Curve3, loc: &glam::DAffine3) -> Curve3 {
    match curve {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: loc.transform_point3(l.origin),
            direction: loc.transform_vector3(l.direction),
        }),
        Curve3::Circle(c) => {
            let center = loc.transform_point3(c.center);
            let normal = loc.transform_vector3(c.normal).normalize_or_zero();
            let x_dir = loc.transform_vector3(c.x_dir).normalize_or_zero();
            let y_dir = loc.transform_vector3(c.y_dir).normalize_or_zero();
            // OCCT-aligned: radius scales by sqrt(scale_in_plane),
            // NOT by normal length.  Use average of x_dir and y_dir
            // transform lengths (in-plane scale factors).
            let sx = loc.transform_vector3(c.x_dir).length().max(1e-12);
            let sy = loc.transform_vector3(c.y_dir).length().max(1e-12);
            let radius = c.radius * (sx * sy).sqrt(); // geometric mean = area scale
            Curve3::Circle(Circle3 {
                center,
                normal,
                x_dir,
                y_dir,
                radius,
            })
        }
        Curve3::BSpline(bs) => Curve3::BSpline(BSplineCurve3 {
            degree: bs.degree,
            knots: bs.knots.clone(),
            control_points: bs
                .control_points
                .iter()
                .map(|&p| loc.transform_point3(p))
                .collect(),
            weights: bs.weights.clone(),
            is_periodic: false,
        }),
        other => other.clone(),
    }
}

/// OCCT-aligned: apply TopLoc_Location transform to a Surface3.
pub fn transform_surface(surface: &Surface3, loc: &glam::DAffine3) -> Surface3 {
    match surface {
        Surface3::Plane(p) => Surface3::Plane(Plane::new(
            loc.transform_point3(p.origin),
            loc.transform_vector3(p.normal).normalize_or_zero(),
        )),
        Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
            origin: loc.transform_point3(c.origin),
            axis: loc.transform_vector3(c.axis).normalize_or_zero(),
            radius: c.radius * loc.transform_vector3(c.axis).length().max(1e-12),
            ref_dir: loc.transform_vector3(c.ref_dir).normalize_or_zero(),
        }),
        Surface3::Sphere(s) => Surface3::Sphere(SphericalSurface {
            center: loc.transform_point3(s.center),
            axis: loc.transform_vector3(s.axis).normalize_or_zero(),
            radius: s.radius * loc.transform_vector3(s.axis).length().max(1e-12),
            ref_dir: loc.transform_vector3(s.ref_dir).normalize_or_zero(),
        }),
        Surface3::BSpline(bs) => Surface3::BSpline(BSplineSurface {
            degree_u: bs.degree_u,
            degree_v: bs.degree_v,
            knots_u: bs.knots_u.clone(),
            knots_v: bs.knots_v.clone(),
            control_points: bs
                .control_points
                .iter()
                .map(|row| row.iter().map(|&p| loc.transform_point3(p)).collect())
                .collect(),
            weights: bs.weights.clone(),
        }),
        other => other.clone(),
    }
}

pub mod eval;
#[cfg(test)]
pub mod tests;
