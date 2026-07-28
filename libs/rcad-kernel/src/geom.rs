use glam::{DVec2, DVec3};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

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
    /// a `gp_Ax3` with the following X direction computation:
    /// 1. Use X (1,0,0) as the reference direction
    /// 2. If |N·X| ≥ 1 - Precision::Angular() (N parallel to X), use Z (0,0,1) instead
    /// 3. u_dir = (ref - N * (N · ref)).normalize()
    /// 4. v_dir = N × u_dir
    pub fn new(origin: DVec3, normal: DVec3) -> Self {
        let normal = normal.normalize_or_zero();
        // OCCT-aligned: gp_Ax3 reference direction selection
        let ref_dir = if normal.x.abs() > 1.0 - 1e-12 {
            DVec3::Z
        } else {
            DVec3::X
        };
        let u_dir = (ref_dir - normal * ref_dir.dot(normal)).normalize_or_zero();
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
}

impl ConicalSurface {
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
        let x_ax = any_perpendicular(axis);
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
}

/// OCCT-aligned: `Geom2d_Conic` intermediate abstract class.
///
/// Groups 2D conic curves (Circle, Ellipse, Hyperbola, Parabola).
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
        // dir_perp forms a right-handed system: axis_dir × normal gives perpendicular direction
        let dir_perp = self.axis_dir.cross(self.normal).normalize();
        self.vertex + (t * t / (2.0 * self.focal_param)) * self.axis_dir + t * dir_perp
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let dir_perp = self.axis_dir.cross(self.normal).normalize();
        let v = (t / self.focal_param) * self.axis_dir + dir_perp;
        v.normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        let dir_perp = self.axis_dir.cross(self.normal).normalize();
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
        let x_ax = any_perpendicular(axis);
        let y_ax = axis.cross(x_ax).normalize();
        let radial = self.radius_at_slant(v);
        let axial = self.axial_from_slant(v);
        self.apex + axial * axis + radial * (u.cos() * x_ax + u.sin() * y_ax)
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let axis = self.axis_dir();
        let x_ax = any_perpendicular(axis);
        let y_ax = axis.cross(x_ax).normalize();
        let radial = u.cos() * x_ax + u.sin() * y_ax;
        let half = self.half_angle_rad;
        (radial * half.cos() - axis * half.sin()).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, f64::INFINITY]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let axis = self.axis_dir();
        let x_ax = any_perpendicular(axis);
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
}

/// OCCT-aligned: `Geom_ElementarySurface` intermediate abstract class.
///
/// Groups analytic surfaces (Plane, Cylinder, Sphere, Cone, Torus, Ellipsoid)
/// with shared access to position, axis, and local frame — corresponding to
/// OCCT's `gp_Ax3` / `gp_Ax2` members stored in each elementary surface.
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
    fn x_axis(&self) -> DVec3 { any_perpendicular(self.axis_dir()) }
    fn y_axis(&self) -> DVec3 { self.axis_dir().cross(any_perpendicular(self.axis_dir())).normalize_or_zero() }
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
fn de_boor_homo(
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
    let mut d: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
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
    d[degree]
}

/// De Boor's algorithm in homogeneous 3D space for 2D rational curves.
/// Returns `[wx, wy, w]` (not divided by w yet).
fn de_boor_homo_2d(
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
    let mut d: Vec<[f64; 3]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights[idx];
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
fn de_boor(degree: usize, knots: &[f64], points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
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
fn de_boor_2d(degree: usize, knots: &[f64], points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
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
    let mut d: Vec<[f64; 3]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights[idx];
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
fn bspline_tangent_analytic(
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
            a_prime.push(s * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
            w_prime.push(DVec3::new(s * (weights[i + 1] - weights[i]), 0.0, 0.0));
        }
    }

    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0f64; m];

    // A'(t): non-rational B-Spline of degree p-1
    let a_prime_t = de_boor(degree - 1, deriv_knots, &a_prime, &unit, t);
    // W'(t): scalar B-Spline of degree p-1 (embedded in .x)
    let w_prime_t = de_boor(degree - 1, deriv_knots, &w_prime, &unit, t).x;

    // W(t) and C(t) from the homogeneous evaluation
    let h = de_boor_homo(degree, knots, points, weights, t);
    let w_t = h[3];
    if w_t.abs() < 1e-15 {
        return DVec3::ZERO;
    }
    let c_t = DVec3::new(h[0] / w_t, h[1] / w_t, h[2] / w_t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

fn bspline_tangent_analytic_2d(
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
            a_prime.push(s * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
            w_prime.push(DVec2::new(s * (weights[i + 1] - weights[i]), 0.0));
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
fn bezier_tangent_analytic(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
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
fn bezier_tangent_analytic_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
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
        de_boor(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        bspline_tangent_analytic(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
        .normalize_or_zero()
    }
    fn derivative_at(&self, t: f64) -> DVec3 {
        bspline_tangent_analytic(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
    fn default_domain(&self) -> [f64; 2] {
        let d = self.degree;
        let n = self.knots.len();
        if n < 2 * d + 2 {
            return [0.0, 1.0];
        }
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
                de_boor_homo(self.degree_u, &self.knots_u, &pts, &wts, u)
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
        de_boor(self.degree_v, &self.knots_v, &v_pts, &v_wts, v)
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
        de_boor_2d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        bspline_tangent_analytic_2d(
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

/// OCCT-aligned: apply TopLoc_Location transform to a Curve3.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle_involute_starts_on_base_circle() {
        let inv = CircleInvolute2d {
            center: DVec2::new(2.0, -1.0),
            base_radius: 3.0,
            start_angle: 0.0,
        };
        let p0 = inv.point_at(0.0);
        assert!((p0.x - 5.0).abs() < 1e-12);
        assert!((p0.y + 1.0).abs() < 1e-12);
    }

    #[test]
    fn archimedean_spiral_point_progresses_radially() {
        let s = ArchimedeanSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.5,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p1 = s.point_at(2.0);
        assert!(
            p1.length() > p0.length(),
            "spiral radius should increase with t"
        );
    }

    #[test]
    fn logarithmic_spiral_grows_exponentially() {
        let s = LogarithmicSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.4,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p1 = s.point_at(2.0);
        assert!(
            p1.length() > p0.length() * 1.5,
            "log spiral should grow faster than linear at this sample"
        );
    }

    #[test]
    fn sine_wave_samples_match_expected_values() {
        let s = SineWave2d {
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p90 = s.point_at(std::f64::consts::FRAC_PI_2);
        assert!((p0.x - 0.0).abs() < 1e-12);
        assert!((p0.y - 0.0).abs() < 1e-12);
        assert!((p90.y - 2.0).abs() < 1e-12);
    }

    #[test]
    fn bspline_curve3_derivative_matches_linear_curve() {
        let curve = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(1.0, 2.0, 3.0), DVec3::new(4.0, 8.0, 15.0)],
            weights: vec![1.0, 1.0],
        };

        let derivative = curve.derivative_at(0.4);

        assert!((derivative - DVec3::new(3.0, 6.0, 12.0)).length() < 1e-12);
    }

    #[test]
    fn bspline_curve2_derivative_matches_linear_curve() {
        let curve = BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::new(1.0, 2.0), DVec2::new(4.0, 8.0)],
            weights: vec![1.0, 1.0],
        };

        let derivative = curve.derivative_at(0.4);

        assert!((derivative - DVec2::new(3.0, 6.0)).length() < 1e-12);
    }

    #[test]
    fn curve2d_sine_wave_variant_dispatches_evaluator() {
        let c = Curve2d::SineWave(SineWave2d {
            amplitude: 1.5,
            frequency: 2.0,
            phase: 0.25,
        });
        let t = 0.3;
        let p = c.point_at(t);
        let expected_y = 1.5 * (2.0 * t + 0.25).sin();
        assert!((p.x - t).abs() < 1e-12);
        assert!((p.y - expected_y).abs() < 1e-12);
    }

    #[test]
    fn sine_wave3_origin_phase_zero_evaluates_at_zero_offset() {
        let c = SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 3.0,
            frequency: 1.0,
            phase: 0.0,
        };
        // At t=0, sin(0)=0 → point should be at origin.
        let p = c.point_at(0.0);
        assert!(
            p.length() < 1e-12,
            "phase-zero at t=0 should be at origin: {p:?}"
        );
        // At t=pi/2, sin(pi/2)=1 → y should equal amplitude.
        let p2 = c.point_at(std::f64::consts::FRAC_PI_2);
        assert!(
            (p2.y - 3.0).abs() < 1e-9,
            "y at t=pi/2 should be amplitude=3: {p2:?}"
        );
    }

    #[test]
    fn curve3_sine_wave_variant_dispatches_evaluator() {
        let c = Curve3::SineWave(SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 1.0,
            frequency: 2.0,
            phase: 0.0,
        });
        let t = 0.5;
        let p = c.point_at(t);
        let expected = DVec3::new(0.5, (2.0_f64 * t).sin(), 0.0);
        assert!((p - expected).length() < 1e-12);
        // Tangent should be non-zero
        let tan = c.tangent_at(t);
        assert!(
            tan.length() > 0.9,
            "tangent should be roughly unit-length: {tan:?}"
        );
    }
}

/// De Casteljau algorithm for rational Bezier curve evaluation in 3D.
/// `t` is in `[0, 1]`.
fn de_casteljau_3d(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
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
fn de_casteljau_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
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

#[cfg(test)]
mod eval_tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn line3_point_at() {
        let l = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        assert!((l.point_at(3.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn circle3_point_at_zero_is_on_circle() {
        // Circle in XY plane, normal = Z
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let p0 = c.point_at(0.0);
        assert!((p0.length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn circle3_full_revolution_closes() {
        let c = Circle3::new(DVec3::new(1.0, 2.0, 3.0), DVec3::Y, 5.0);
        let p0 = c.point_at(0.0);
        let p2pi = c.point_at(2.0 * PI);
        assert!((p0 - p2pi).length() < 1e-10);
    }

    #[test]
    fn circle3_quarter_turn() {
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let p0 = c.point_at(0.0);
        let p90 = c.point_at(FRAC_PI_2);
        // 90° rotation: p0 and p90 should be perpendicular from center
        assert!((p0.dot(p90)).abs() < 1e-10);
        assert!((p90.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_surface_north_pole() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
            ref_dir: any_perpendicular(DVec3::Y),
        };
        // v=0 is north pole regardless of u
        let p = s.point_at(0.0, 0.0);
        // Should be at (0, 3, 0)
        assert!((p - DVec3::new(0.0, 3.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn sphere_surface_point_on_sphere() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Y),
        };
        for u in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            for v in [0.1, 0.5, 1.0, PI / 2.0, PI - 0.1] {
                let p = s.point_at(u, v);
                assert!(
                    (p.length() - 2.0).abs() < 1e-9,
                    "u={u} v={v} |p|={}",
                    p.length()
                );
            }
        }
    }

    #[test]
    fn cylinder_surface_point_on_cylinder() {
        let c = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        for u in [0.0, 1.0, PI, 2.0 * PI - 0.1] {
            let p = c.point_at(u, 0.0);
            let radial = DVec3::new(p.x, 0.0, p.z).length();
            assert!((radial - 3.0).abs() < 1e-9, "u={u} radial={radial}");
        }
    }

    #[test]
    fn bspline_degree1_linear_interpolation() {
        // Degree-1 BSpline with 2 control points = straight line
        let c = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };
        let p0 = c.point_at(0.0);
        let p1 = c.point_at(1.0);
        let pmid = c.point_at(0.5);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);
        assert!((p1 - DVec3::X).length() < 1e-10);
        assert!((pmid - DVec3::new(0.5, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn bspline_degree2_quadratic() {
        // Degree-2 quadratic arc through 3 control points
        let c = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::new(0.5, 1.0, 0.0), DVec3::X],
            weights: vec![1.0, 1.0, 1.0],
        };
        let p0 = c.point_at(0.0);
        let p1 = c.point_at(1.0);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);
        assert!((p1 - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn torus_surface_point_on_torus() {
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        for u in [0.0, PI / 2.0, PI] {
            for v in [0.0, PI / 2.0, PI] {
                let p = t.point_at(u, v);
                // Distance from the tube center circle should be minor_radius
                let x_ax = any_perpendicular(DVec3::Y);
                let y_ax = DVec3::Y.cross(x_ax).normalize();
                let tube_center = t.center + t.major_radius * (u.cos() * x_ax + u.sin() * y_ax);
                assert!((p - tube_center).length() - 1.0 < 1e-9, "u={u} v={v}");
            }
        }
    }

    #[test]
    fn ellipsoid_surface_satisfies_implicit_equation() {
        let s = EllipsoidalSurface {
            center: DVec3::new(1.0, -2.0, 0.5),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 4.0,
            radius_y: 2.0,
            radius_z: 1.5,
        };
        let p = s.point_at(0.7, 1.2) - s.center;
        let value =
            (p.x / s.radius_x).powi(2) + (p.y / s.radius_y).powi(2) + (p.z / s.radius_z).powi(2);
        assert!(
            (value - 1.0).abs() < 1e-9,
            "implicit value should be 1, got {value}"
        );
    }

    #[test]
    fn ellipsoid_surface_normal_matches_gradient_direction() {
        let s = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };
        let u = 0.9;
        let v = 1.1;
        let p = s.point_at(u, v);
        let expected = DVec3::new(
            p.x / (s.radius_x * s.radius_x),
            p.y / (s.radius_y * s.radius_y),
            p.z / (s.radius_z * s.radius_z),
        )
        .normalize();
        let n = s.normal_at(u, v);
        assert!(
            (n - expected).length() < 1e-9,
            "n={n:?} expected={expected:?}"
        );
    }

    #[test]
    fn helicoid_surface_advances_by_pitch_per_turn() {
        let s = HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 6.0,
        };
        let p0 = s.point_at(0.0, 2.0);
        let p1 = s.point_at(2.0 * PI, 2.0);
        let delta = p1 - p0;
        assert!(
            (delta - DVec3::new(0.0, 0.0, 6.0)).length() < 1e-9,
            "delta={delta:?}"
        );
    }

    #[test]
    fn helicoid_surface_normal_is_perpendicular_to_parametric_tangents() {
        let s = HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 4.0,
        };
        let u = 0.6;
        let v = 1.75;
        let n = s.normal_at(u, v);
        let eps = 1e-6;
        let du = (s.point_at(u + eps, v) - s.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (s.point_at(u, v + eps) - s.point_at(u, v - eps)) / (2.0 * eps);
        assert!(
            n.dot(du).abs() < 1e-6,
            "n·du={} should be near 0",
            n.dot(du)
        );
        assert!(
            n.dot(dv).abs() < 1e-6,
            "n·dv={} should be near 0",
            n.dot(dv)
        );
        assert!(n.length() > 0.99, "normal should be unit-length: {n:?}");
    }

    #[test]
    fn pipe_surface_with_line_spine_matches_cylindrical_section() {
        let surface = PipeSurface {
            spine: Box::new(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::Z,
            })),
            ref_dir: DVec3::X,
            radius: 2.0,
        };

        assert!((surface.point_at(0.0, 0.0) - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);
        assert!((surface.point_at(PI * 0.5, 0.5) - DVec3::new(0.0, 2.0, 0.5)).length() < 1e-9);
        assert!((surface.default_domain()[0] - 0.0).abs() < 1e-12);
        assert!((surface.default_domain()[1] - 2.0 * PI).abs() < 1e-12);
    }

    #[test]
    fn tri_bezier_surface_hits_triangle_corners() {
        let surface = TriBezierSurface {
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0]],
        };
        assert!((surface.point_at(0.0, 0.0) - DVec3::new(0.0, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(1.0, 0.0) - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.0, 1.0) - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn tri_bezier_surface_dispatches_through_surface3() {
        let surface = Surface3::TriBezier(TriBezierSurface {
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0]],
        });
        let p = surface.point_at(0.25, 0.5);
        assert!(p.x >= -1e-12 && p.y >= -1e-12);
        assert!(surface.normal_at(0.2, 0.2).length() > 0.99);
    }

    #[test]
    fn ruled_surface_interpolates_between_curves() {
        let surface = RuledSurface {
            start: Box::new(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::X,
            })),
            end: Box::new(Curve3::Line(Line3 {
                origin: DVec3::Y,
                direction: DVec3::X,
            })),
        };
        assert!((surface.point_at(0.25, 0.0) - DVec3::new(0.25, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.25, 1.0) - DVec3::new(0.25, 1.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.25, 0.5) - DVec3::new(0.25, 0.5, 0.0)).length() < 1e-12);
        assert!(surface.normal_at(0.25, 0.5).length() > 0.99);
    }

    #[test]
    fn coons_surface_interpolates_all_four_boundaries() {
        let surface = CoonsSurface {
            south: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::X,
            })),
            north: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 1.0, 1.0),
                direction: DVec3::X,
            })),
            west: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
            east: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(1.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
        };

        assert!((surface.point_at(0.3, 0.0) - DVec3::new(0.3, 0.0, 0.0)).length() < 1e-9);
        assert!((surface.point_at(0.3, 1.0) - DVec3::new(0.3, 1.0, 1.0)).length() < 1e-9);
        assert!((surface.point_at(0.0, 0.4) - DVec3::new(0.0, 0.4, 0.4)).length() < 1e-9);
        assert!((surface.point_at(1.0, 0.4) - DVec3::new(1.0, 0.4, 0.4)).length() < 1e-9);
        assert!((surface.point_at(0.5, 0.5) - DVec3::new(0.5, 0.5, 0.5)).length() < 1e-9);
    }

    #[test]
    fn conical_surface_uses_slant_distance_from_reference_circle() {
        let surface = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 30.0_f64.to_radians(),
        };

        let p0 = surface.point_at(0.0, 0.0);
        assert!(p0.dot(surface.axis_dir()).abs() < 1e-9);
        assert!((p0.length() - 2.0).abs() < 1e-9);

        let slant = 4.0;
        let p1 = surface.point_at(0.0, slant);
        assert!((p1.z - slant * surface.half_angle_rad.cos()).abs() < 1e-9);
        let radial = p1 - surface.axis_dir() * p1.dot(surface.axis_dir());
        assert!((radial.length() - (2.0 + slant * surface.half_angle_rad.sin())).abs() < 1e-9);
    }

    #[test]
    fn conical_surface_derives_true_apex_from_reference_circle() {
        let surface = ConicalSurface {
            apex: DVec3::new(0.0, 0.0, 5.0),
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 45.0_f64.to_radians(),
        };

        assert!((surface.apex_point() - DVec3::new(0.0, 0.0, 3.0)).length() < 1e-9);
    }

    // --- Analytic derivative tests ---

    /// Quadratic Bezier: P0=(0,0,0), P1=(0.5,1,0), P2=(1,0,0), unit weights.
    /// Analytic tangent at t=0 should be (0.5,1,0).normalize() = (1,2,0)/√5.
    #[test]
    fn bezier_tangent_at_endpoint_analytic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(0.5, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let c = BezierCurve3 {
            control_points: pts,
            weights: wts,
        };
        let tan = c.tangent_at(0.0);
        let expected = DVec3::new(1.0, 2.0, 0.0).normalize();
        assert!(
            (tan - expected).length() < 1e-10,
            "tan={tan:?} expected={expected:?}"
        );
    }

    /// Quadratic Bezier tangent at t=1 should be (1,-2,0)/√5.
    #[test]
    fn bezier_tangent_at_end_analytic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(0.5, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let c = BezierCurve3 {
            control_points: pts,
            weights: wts,
        };
        let tan = c.tangent_at(1.0);
        let expected = DVec3::new(1.0, -2.0, 0.0).normalize();
        assert!(
            (tan - expected).length() < 1e-10,
            "tan={tan:?} expected={expected:?}"
        );
    }

    /// Degree-1 B-Spline (polyline): tangent should be constant along each segment.
    #[test]
    fn bspline_degree1_tangent_is_segment_direction() {
        // Two-segment polyline: (0,0,0) -> (1,0,0) -> (1,1,0)
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.5, 1.0, 1.0];
        let c = BSplineCurve3 {
            degree: 1,
            knots,
            control_points: pts,
            weights: wts,
        };
        let tan0 = c.tangent_at(0.1);
        assert!(
            (tan0 - DVec3::X).length() < 1e-10,
            "first segment should be +X, got {tan0:?}"
        );
        let tan1 = c.tangent_at(0.9);
        assert!(
            (tan1 - DVec3::Y).length() < 1e-10,
            "second segment should be +Y, got {tan1:?}"
        );
    }

    /// Degree-2 B-Spline circle arc: tangent should be perpendicular to radius.
    #[test]
    fn bspline_circle_tangent_perpendicular_to_radius() {
        // Use circle_to_bspline to get an exact NURBS circle, then check tangents.
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let c = crate::nurbs_convert::circle_to_bspline(&circle);
        for &t in &[0.0, 0.5, 1.0, 1.5, 2.0] {
            let pt = c.point_at(t);
            let tan = c.tangent_at(t);
            // Tangent must be perpendicular to the radius vector
            let dot = pt.normalize_or_zero().dot(tan);
            assert!(
                dot.abs() < 1e-8,
                "t={t}: radius*tangent={dot} (should be 0)"
            );
            // Tangent must be a unit vector
            assert!(
                (tan.length() - 1.0).abs() < 1e-10,
                "t={t}: |tan|={}",
                tan.length()
            );
        }
    }

    // =========================================================================
    // OCCT-aligned TKG3d / TKG2d evaluation tests
    // =========================================================================
    //
    // These test point_at (D0) and tangent_at (D1) for each curve/surface type,
    // matching patterns in OCCT's TKG3d/GTests/ and TKG2d/GTests/.

    // ── 3D Curve evaluation (Geom_CurveEval_Test.cxx pattern) ────────────

    #[test]
    fn line_eval_d0_d1() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let p = line.point_at(5.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-12);
        let t = line.tangent_at(5.0);
        assert!((t - DVec3::X).length() < 1e-12);
    }

    #[test]
    fn line_eval_d2_zero_second_derivative() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        // For a line, the first derivative (tangent) is constant; second derivative is zero.
        // The curve is linear: P(t) = origin + t * direction
        // The second derivative: d²P/dt² = 0
        // Using the derivative_at method:
        let d1 = line.derivative_at(0.0);
        let d2 = (line.derivative_at(1e-4) - line.derivative_at(-1e-4)) / (2.0 * 1e-4);
        assert!((d1 - DVec3::X).length() < 1e-10);
        assert!(
            d2.length() < 1e-10,
            "Line second derivative should be 0, got {d2:?}"
        );
    }

    #[test]
    fn circle_eval_d0_d1() {
        // Circle3::new(ZERO, Z, 5.0): OCCT gp_Ax2 gives x_dir=X, y_dir=Z×X=Y
        // P(0) = 5*X = (5,0,0), tangent = Y
        // P(PI/2) = 5*Y = (0,5,0), tangent = -X
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let p0 = circle.point_at(0.0);
        assert!((p0 - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-10);
        let p_half = circle.point_at(std::f64::consts::PI / 2.0);
        assert!((p_half - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10);
        let p_pi = circle.point_at(std::f64::consts::PI);
        assert!((p_pi - DVec3::new(-5.0, 0.0, 0.0)).length() < 1e-10);
        // Tangent at 0 should be (0, 1, 0)
        let t0 = circle.tangent_at(0.0);
        assert!((t0 - DVec3::Y).length() < 1e-10);
        // Tangent at PI/2 should be (-1, 0, 0)
        let t_half = circle.tangent_at(std::f64::consts::PI / 2.0);
        assert!((t_half - DVec3::new(-1.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn circle_transform_copy() {
        // Circle3::new((1,2,3), Z, 4.0): OCCT gp_Ax2 gives x_dir=X, y_dir=Y
        // P(0) = (1,2,3) + 4*X = (5,2,3)
        let circle = Circle3::new(DVec3::new(1.0, 2.0, 3.0), DVec3::Z, 4.0);
        assert!((circle.point_at(0.0) - DVec3::new(5.0, 2.0, 3.0)).length() < 1e-10);
    }

    #[test]
    fn bspline_eval_d0_d1_d2_consistency() {
        // Degree-2 BSpline through 4 points
        let c = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.3, 0.7, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(2.0, 3.0, 0.0),
                DVec3::new(5.0, 3.0, 0.0),
                DVec3::new(7.0, 0.0, 0.0),
                DVec3::new(10.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 5],
        };
        // Point at t=0 should be first control point
        assert!((c.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // Point at t=1 should be last control point
        assert!((c.point_at(1.0) - DVec3::new(10.0, 0.0, 0.0)).length() < 1e-10);
        // Range check: some midpoint
        let _pmid = c.point_at(0.5);
        assert!(_pmid.x >= 0.0 && _pmid.x <= 10.0);
    }

    #[test]
    fn bezier_eval_d0_d1() {
        // Cubic Bezier through 4 control points
        let c = BezierCurve3 {
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 3.0, 0.0),
                DVec3::new(3.0, 3.0, 0.0),
                DVec3::new(4.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        };
        // t=0 -> first pole
        assert!((c.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // t=1 -> last pole
        assert!((c.point_at(1.0) - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-10);
        // Tangent at t=0: direction from P0 to P1
        let t0 = c.tangent_at(0.0);
        assert!((t0 - DVec3::new(1.0, 3.0, 0.0).normalize()).length() < 1e-10);
    }

    // ── 3D Surface evaluation (Geom_SurfaceEval_Test.cxx pattern) ───────

    #[test]
    fn plane_eval_d0_d1() {
        // Plane with normal Z: OCCT gp_Ax3 gives u_dir=X, v_dir=Y.
        // P(u,v) = u*X + v*Y. So P(2,3) = (2, 3, 0)
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let p = plane.point_at(2.0, 3.0);
        // Should be in the plane (Z=0)
        assert!((p.z).abs() < 1e-10);
        // Distance from origin should be sqrt(4+9) = sqrt(13)
        assert!((p.length() - 13.0f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn cylinder_eval_d0() {
        // Cylinder with axis Z, any_perpendicular(Z) = Y, y_ax = Z.cross(Y) = -X
        // Cylinder with axis Z, ref_dir=X:
        // x_ax = X, y_ax = Z×X = Y
        // P(u,v) = R*(cos(u)*X + sin(u)*Y) + v*Z
        // P(0,0) = 3*X = (3, 0, 0)
        // P(PI/2, 5) = 3*Y + 5*Z = (0, 3, 5)
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        let p0 = cyl.point_at(0.0, 0.0);
        assert!((p0 - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
        let p1 = cyl.point_at(std::f64::consts::PI / 2.0, 5.0);
        assert!((p1 - DVec3::new(0.0, 3.0, 5.0)).length() < 1e-10);
    }

    #[test]
    fn sphere_eval_d0_full_sphere() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: DVec3::X,
        };
        for u in [0.0, 1.0, 2.0, 4.0, 6.0] {
            for v in [-1.0, -0.5, 0.0, 0.5, 1.0] {
                let p = s.point_at(u, v);
                assert!((p.length() - 2.0).abs() < 1e-9, "u={u} v={v}");
            }
        }
    }

    #[test]
    fn cone_eval_d0() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 45.0_f64.to_radians(),
        };
        // At V=0, radius=2
        let p0 = cone.point_at(0.0, 0.0);
        assert!((p0.x - 2.0).abs() < 1e-9 || (p0.y - 2.0).abs() < 1e-9);
        assert!((p0.z).abs() < 1e-9);
    }

    #[test]
    fn torus_eval_d0() {
        use std::f64::consts::PI;
        // Torus with axis Z: OCCT-aligned x_ax = X, y_ax = Y
        // P(u,v) = (R + r*cos(v))*(cos(u)*X + sin(u)*Y) + r*sin(v)*Z
        // P(0,0) = (5+1)*X = (6, 0, 0)
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let p = t.point_at(0.0, 0.0);
        assert!((p - DVec3::new(6.0, 0.0, 0.0)).length() < 1e-9);
        // At u=0, v=PI: P = (5-1)*X = (4, 0, 0)
        let p2 = t.point_at(0.0, PI);
        assert!((p2 - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-9);
    }

    #[test]
    fn offset_surface_eval_d0() {
        // Use a Sphere as basis — normal is well-defined
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 1.0);
        let off = OffsetSurface {
            basis: Box::new(Surface3::Sphere(sphere)),
            offset_distance: 0.5,
        };
        // Offset sphere: radius = 1 + 0.5 = 1.5
        let p = off.point_at(0.0, 0.0);
        assert!((p.length() - 1.5).abs() < 1e-9);
    }

    // ── 2D Curve evaluation (Geom2d_CurveEval_Test.cxx pattern) ─────────

    #[test]
    fn line2d_eval_d0() {
        let l = Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        };
        let p = l.point_at(3.0);
        assert!((p - DVec2::new(3.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn circle2d_eval_d0() {
        use std::f64::consts::PI;
        let c = Circle2d::new(DVec2::ZERO, 5.0);
        let p0 = c.point_at(0.0);
        assert!((p0 - DVec2::new(5.0, 0.0)).length() < 1e-12);
        let p_half = c.point_at(PI / 2.0);
        assert!((p_half - DVec2::new(0.0, 5.0)).length() < 1e-12);
        let p_pi = c.point_at(PI);
        assert!((p_pi - DVec2::new(-5.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn circle2d_revolved() {
        let c = Circle2d::new(DVec2::ZERO, 5.0);
        let p0 = c.point_at(0.0);
        let p2pi = c.point_at(2.0 * std::f64::consts::PI);
        assert!((p0 - p2pi).length() < 1e-12);
    }

    #[test]
    fn ellipse2d_eval_d0() {
        use std::f64::consts::PI;
        let e = Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        let p0 = e.point_at(0.0);
        assert!((p0 - DVec2::new(10.0, 0.0)).length() < 1e-12);
        let p_half = e.point_at(PI / 2.0);
        assert!((p_half - DVec2::new(0.0, 5.0)).length() < 1e-12);
    }

    #[test]
    fn parabola2d_eval_d0() {
        let p = Parabola2d {
            origin: DVec2::ZERO,
            axis_dir: DVec2::X,
            focal_param: 4.0,
        };
        // P(t) = (t²/(2p), t) = (t²/8, t)
        let p0 = p.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-12);
        let p4 = p.point_at(4.0);
        assert!((p4 - DVec2::new(2.0, 4.0)).length() < 1e-10);
    }

    #[test]
    fn hyperbola2d_eval_d0() {
        let h = Hyperbola2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };
        let p0 = h.point_at(0.0);
        // X = a*cosh(0) = 3, Y = b*sinh(0) = 0
        assert!((p0 - DVec2::new(3.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn bspline2_eval_d0() {
        let c = BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::ZERO, DVec2::X],
            weights: vec![1.0, 1.0],
        };
        assert!((c.point_at(0.0) - DVec2::ZERO).length() < 1e-12);
        assert!((c.point_at(1.0) - DVec2::X).length() < 1e-12);
        assert!((c.point_at(0.5) - DVec2::new(0.5, 0.0)).length() < 1e-12);
    }

    #[test]
    fn bezier2_eval_d0() {
        let c = BezierCurve2 {
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(0.5, 1.0),
                DVec2::new(1.0, 0.0),
            ],
            weights: vec![1.0; 3],
        };
        assert!((c.point_at(0.0) - DVec2::ZERO).length() < 1e-12);
        assert!((c.point_at(1.0) - DVec2::X).length() < 1e-12);
    }

    #[test]
    fn trimmed_curve2_eval() {
        let inner = Circle2d::new(DVec2::ZERO, 5.0);
        let tc = TrimmedCurve2 {
            curve: Box::new(Curve2d::Circle(inner)),
            t_min: 0.0,
            t_max: std::f64::consts::PI,
        };
        let p0 = tc.point_at(0.0);
        assert!((p0 - DVec2::new(5.0, 0.0)).length() < 1e-12);
        let p_half = tc.point_at(std::f64::consts::PI / 2.0);
        assert!((p_half - DVec2::new(0.0, 5.0)).length() < 1e-12);
        // Out of range → clamped
        let p_out = tc.point_at(10.0);
        assert!((p_out - DVec2::new(-5.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn offset_curve2d_eval() {
        // Circle2d r=5, offset_distance uses right-hand normal.
        // P(0) = (5,0) + 1*(right_normal at 0) = (5,0) + 1*(1,0) ≈ (6,0)
        let basis = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let off = OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 1.0,
        };
        let p0 = off.point_at(0.0);
        // Offset point should be outside the original circle
        assert!(p0.length() > 5.0);
        assert!((p0.length() - 6.0).abs() < 0.1);
    }

    // ── Special 2D curve evaluation tests ───────────────────────────────

    #[test]
    fn circle_involute2d_eval() {
        let inv = CircleInvolute2d {
            center: DVec2::ZERO,
            base_radius: 3.0,
            start_angle: 0.0,
        };
        // At t=0: P = center + r*(cos(0)+0*sin(0), sin(0)-0*cos(0)) = center + (r, 0)
        let p0 = inv.point_at(0.0);
        assert!((p0 - DVec2::new(3.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn sine_wave2d_eval() {
        let w = SineWave2d {
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        };
        let p0 = w.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-12);
        let p_pi = w.point_at(std::f64::consts::PI / 2.0);
        assert!((p_pi - DVec2::new(std::f64::consts::PI / 2.0, 2.0)).length() < 1e-10);
    }

    #[test]
    fn archimedean_spiral2d_eval() {
        let s = ArchimedeanSpiral2d {
            center: DVec2::ZERO,
            a: 0.0,
            b: 1.0,
            start_angle: 0.0,
        };
        // r(t) = 0 + 1*t, theta(t) = t
        let p0 = s.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-12);
        let p1 = s.point_at(2.0);
        assert!((p1.length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn logarithmic_spiral2d_eval() {
        let s = LogarithmicSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.5,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        assert!((p0 - DVec2::new(1.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn offset_curve3_eval() {
        // Circle3 offset along Z — FD tangent makes this approximate
        let basis = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let off = OffsetCurve3 {
            basis: Box::new(basis),
            offset_distance: 2.0,
            offset_dir: DVec3::Z,
        };
        let p0 = off.point_at(0.0);
        // With offset, should differ from original circle (radius 5)
        assert!((p0.length() - 5.0).abs() > 0.5);
        // But not be wildly different
        assert!(p0.length() < 10.0);
    }

    // ── Reverse parameter tests (Geom2d_*_ReversedParameter pattern) ────

    #[test]
    fn circle_reverse_eval() {
        // A reversed circle should evaluate in the opposite direction
        // OCCT: ReversedParameter(t) = 2*PI - t
        // Since rcad doesn't have an explicit reverse flag, verify that
        // the parameterization wraps around correctly.
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let p_fwd = c.point_at(0.3);
        let p_rev = c.point_at(2.0 * std::f64::consts::PI - 0.3);
        // These are different points (one advances forward, one backward)
        assert!((p_fwd - p_rev).length() > 0.5);
    }

    // =========================================================================
    // OCCT-aligned comprehensive 3D curve evaluation tests
    // (matching TKG3d/GTests Geom_Line/Circle/Ellipse/BSpline/Bezier patterns)
    // =========================================================================

    // ── Line ────────────────────────────────────────────────────────────

    #[test]
    fn line3_eval_at_multiple_points() {
        let line = Line3 {
            origin: DVec3::new(1.0, 2.0, 3.0),
            direction: DVec3::new(0.0, 1.0, 0.0),
        };
        assert!((line.point_at(0.0) - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-12);
        assert!((line.point_at(5.0) - DVec3::new(1.0, 7.0, 3.0)).length() < 1e-12);
        assert!((line.point_at(-3.0) - DVec3::new(1.0, -1.0, 3.0)).length() < 1e-12);
    }

    #[test]
    fn line3_constant_tangent_and_derivative() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::new(1.0, 1.0, 1.0).normalize(),
        };
        let d = DVec3::new(1.0, 1.0, 1.0).normalize();
        for &t in &[-10.0, -1.0, 0.0, 1.0, 10.0] {
            assert!((line.tangent_at(t) - d).length() < 1e-12);
            assert!((line.derivative_at(t) - d).length() < 1e-12);
        }
    }

    #[test]
    fn line3_default_domain_infinite() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let [t0, t1] = line.default_domain();
        assert!(t0.is_infinite() && t0.is_sign_negative());
        assert!(t1.is_infinite() && t1.is_sign_positive());
    }

    // ── Circle ─────────────────────────────────────────────────────────

    #[test]
    fn circle3_eval_four_quadrants() {
        // Use explicit x_dir/y_dir for predictable orientation
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            radius: 5.0,
        };
        // P(0) = (5,0,0), P(PI/2) = (0,5,0), P(PI) = (-5,0,0), P(3PI/2) = (0,-5,0)
        assert!((c.point_at(0.0) - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-10);
        assert!(
            (c.point_at(std::f64::consts::PI / 2.0) - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10
        );
        assert!((c.point_at(std::f64::consts::PI) - DVec3::new(-5.0, 0.0, 0.0)).length() < 1e-10);
        assert!(
            (c.point_at(3.0 * std::f64::consts::PI / 2.0) - DVec3::new(0.0, -5.0, 0.0)).length()
                < 1e-10
        );
    }

    #[test]
    fn circle3_tangent_at_quadrants() {
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            radius: 5.0,
        };
        use std::f64::consts::PI;
        // tangent = (-R*sin(t)*X + R*cos(t)*Y).normalize()
        assert!((c.tangent_at(0.0) - DVec3::Y).length() < 1e-10);
        assert!((c.tangent_at(PI / 2.0) + DVec3::X).length() < 1e-10);
        assert!((c.tangent_at(PI) + DVec3::Y).length() < 1e-10);
        assert!((c.tangent_at(3.0 * PI / 2.0) - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn circle3_derivative_nonzero() {
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            radius: 5.0,
        };
        // derivative = R * (-sin(t)*X + cos(t)*Y), always non-zero for R>0
        for &t in &[0.0, 0.5, 1.0, 2.0, 4.0, 6.0] {
            let d = c.derivative_at(t);
            assert!(d.length() > 0.0);
            assert!(
                (d.length() - 5.0).abs() < 1e-10,
                "t={} |d|={}",
                t,
                d.length()
            );
        }
    }

    #[test]
    fn circle3_default_domain() {
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let [t0, t1] = c.default_domain();
        assert!((t0 - 0.0).abs() < 1e-12);
        assert!((t1 - 2.0 * std::f64::consts::PI).abs() < 1e-12);
    }

    // ── Ellipse ────────────────────────────────────────────────────────

    #[test]
    fn ellipse3_eval_vertices() {
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        use std::f64::consts::PI;
        // Major vertices: t=0 → (10,0,0), t=PI → (-10,0,0)
        assert!((e.point_at(0.0) - DVec3::new(10.0, 0.0, 0.0)).length() < 1e-10);
        assert!((e.point_at(PI) - DVec3::new(-10.0, 0.0, 0.0)).length() < 1e-10);
        // Minor vertices: t=PI/2 → (0,5,0), t=3PI/2 → (0,-5,0)
        assert!((e.point_at(PI / 2.0) - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10);
        assert!((e.point_at(3.0 * PI / 2.0) - DVec3::new(0.0, -5.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn ellipse3_tangent_at_major_vertex() {
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        // Tangent at major vertex t=0: direction = (0, 5, 0) = Y
        // (derivative: -a*sin(0)*X + b*cos(0)*Y = 5*Y)
        let t0 = e.tangent_at(0.0);
        assert!((t0 - DVec3::Y).length() < 1e-10);
    }

    #[test]
    fn ellipse3_derivative_at_vertices() {
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        // derivative at t=0: a*(-sin(0))*X + b*cos(0)*Y = b*Y = 5*Y
        let d0 = e.derivative_at(0.0);
        assert!((d0 - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10);
        // derivative at t=PI/2: a*(-sin(PI/2))*X + b*cos(PI/2)*Y = -a*X = -10*X
        let d_half = e.derivative_at(std::f64::consts::PI / 2.0);
        assert!((d_half - DVec3::new(-10.0, 0.0, 0.0)).length() < 1e-10);
    }

    // ── BSpline ────────────────────────────────────────────────────────

    #[test]
    fn bspline3_eval_at_knots() {
        let c = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.3, 0.7, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(3.0, 5.0, 0.0),
                DVec3::new(6.0, 5.0, 0.0),
                DVec3::new(9.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        };
        // Endpoints
        assert!((c.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        assert!((c.point_at(1.0) - DVec3::new(9.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn bspline3_degree1_is_line() {
        let c = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(1.0, 2.0, 3.0), DVec3::new(4.0, 5.0, 6.0)],
            weights: vec![1.0, 1.0],
        };
        assert!((c.point_at(0.0) - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-10);
        assert!((c.point_at(1.0) - DVec3::new(4.0, 5.0, 6.0)).length() < 1e-10);
        assert!((c.point_at(0.5) - DVec3::new(2.5, 3.5, 4.5)).length() < 1e-10);
    }

    // ── Bezier ─────────────────────────────────────────────────────────

    #[test]
    fn bezier3_linear_tangent_constant() {
        // Degree-1 Bezier (line) has constant tangent
        let c = BezierCurve3 {
            control_points: vec![DVec3::ZERO, DVec3::new(2.0, 4.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let t0 = c.tangent_at(0.0);
        let t1 = c.tangent_at(0.5);
        let t_end = c.tangent_at(1.0);
        assert!((t0 - t1).length() < 1e-10);
        assert!((t0 - t_end).length() < 1e-10);
    }

    #[test]
    fn bezier3_rational_weight_effect() {
        // Rational Bezier with center weight > 1 pulls curve toward control point
        let non_rational = BezierCurve3 {
            control_points: vec![
                DVec3::new(-1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let rational = BezierCurve3 {
            control_points: vec![
                DVec3::new(-1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 5.0, 1.0], // heavy center weight pulls up
        };
        let p_nr = non_rational.point_at(0.5);
        let p_r = rational.point_at(0.5);
        // Rational with heavy center weight should be higher (more Y)
        assert!(p_r.y > p_nr.y);
    }

    // ── Parabola ───────────────────────────────────────────────────────

    #[test]
    fn parabola3_eval_and_derivative() {
        use std::f64::consts::PI;
        // Parabola: dir_perp = axis_dir.cross(normal) = X.cross(Z) = -Y
        // P(t) = (t²/(2p)) * X + t * (-Y) = (t²/8, -t, 0)
        let p = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 4.0,
        };
        // P(0) = (0,0,0)
        assert!((p.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // P(4) = (16/8, -4, 0) = (2, -4, 0)
        assert!((p.point_at(4.0) - DVec3::new(2.0, -4.0, 0.0)).length() < 1e-10);
        // derivative = (t/4, -1, 0): at t=0 → (0, -1, 0), at t=4 → (1, -1, 0)
        let d0 = p.derivative_at(0.0);
        assert!((d0 - DVec3::new(0.0, -1.0, 0.0)).length() < 1e-10);
        let d4 = p.derivative_at(4.0);
        assert!((d4 - DVec3::new(1.0, -1.0, 0.0)).length() < 1e-10);
    }

    // ── Hyperbola ──────────────────────────────────────────────────────

    #[test]
    fn hyperbola3_eval_and_derivative() {
        use std::f64::consts::PI;
        let h = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };
        // P(0) = (a*cosh(0), 0, 0) = (3, 0, 0)  (since minor_dir = normal×major_dir = Z×X = Y)
        // Actually: P(t) = center + a*cosh(t)*X + b*sinh(t)*Y
        // P(0) = (3*1, 2*0, 0) = (3, 0, 0)
        assert!((h.point_at(0.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
        // P(1) = (3*cosh(1), 2*sinh(1), 0)
        let p1 = h.point_at(1.0);
        assert!((p1.x - 3.0 * 1.0f64.cosh()).abs() < 1e-10);
        assert!((p1.y - 2.0 * 1.0f64.sinh()).abs() < 1e-10);
        // derivative = (a*sinh(t), b*cosh(t), 0)
        // at t=0: (0, 2, 0)
        let d0 = h.derivative_at(0.0);
        assert!((d0 - DVec3::new(0.0, 2.0, 0.0)).length() < 1e-10);
    }

    // ── Helix ──────────────────────────────────────────────────────────

    #[test]
    fn helix3_full_turn_pitch_advance() {
        // Helix with pitch 2: after one full turn, Z advances by 2
        let h = CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 3.0,
            pitch: 2.0,
        };
        use std::f64::consts::{PI, TAU};
        // x_axis = ref_dir - axis * dot = X - 0 = X
        // y_axis = Z.cross(X) = Y
        // P(0) = (3, 0, 0), P(TAU) = (3, 0, pitch) = (3, 0, 2)
        assert!((h.point_at(0.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
        let p_full = h.point_at(TAU);
        assert!((p_full - DVec3::new(3.0, 0.0, 2.0)).length() < 1e-10);
    }

    #[test]
    fn helix3_half_turn_opposite() {
        let h = CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 3.0,
            pitch: 2.0,
        };
        use std::f64::consts::PI;
        // P(PI) = (-3, 0, 1)
        let p_half = h.point_at(PI);
        assert!((p_half - DVec3::new(-3.0, 0.0, 1.0)).length() < 1e-10);
    }

    // ── SineWave ───────────────────────────────────────────────────────

    #[test]
    fn sine_wave3_eval_and_derivative() {
        let w = SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        };
        use std::f64::consts::PI;
        // P(0) = (0, 0, 0)
        assert!((w.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // P(PI/2) = (PI/2, 2, 0)
        let p = w.point_at(PI / 2.0);
        assert!((p - DVec3::new(PI / 2.0, 2.0, 0.0)).length() < 1e-10);
        // derivative = X + 2*cos(t)*Y, at t=0: X + 2*Y
        let d0 = w.derivative_at(0.0);
        assert!((d0 - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-10);
    }

    // ── OffsetCurve3 ───────────────────────────────────────────────────

    #[test]
    fn offset_curve3_line_offset() {
        // Line along X, offset along Z: tangent = X, perp = X×Z = -Y
        // The offset displaces in the -Y direction (perpendicular to both tangent and offset_dir)
        // FD tangent gives approximate direction, so just check the point differs from the line
        let basis = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let off = OffsetCurve3 {
            basis: Box::new(basis),
            offset_distance: 2.0,
            offset_dir: DVec3::Z,
        };
        let pt = off.point_at(5.0);
        // Should differ from the base line point (5,0,0)
        assert!((pt - DVec3::new(5.0, 0.0, 0.0)).length() > 1.0);
        // The Z coordinate should be near 0 (offset in XY plane, not Z)
        assert!(pt.z.abs() < 0.1);
    }

    // =========================================================================
    // OCCT-aligned comprehensive surface evaluation tests
    // (matching TKG3d/GTests Geom_Plane/Cylinder/Sphere/Cone/Torus patterns)
    // =========================================================================

    // ── Plane ───────────────────────────────────────────────────────────

    #[test]
    fn plane_derivatives_constant() {
        let p = Plane::new(DVec3::ZERO, DVec3::Z);
        // OCCT-aligned: gp_Ax3(gp_Pnt, gp_Dir) with normal=Z gives u_dir=X, v_dir=Y
        let (pt, dpu, dpv) = p.derivatives(2.0, 3.0);
        assert!((dpu - DVec3::X).length() < 1e-10);
        assert!((dpv - DVec3::Y).length() < 1e-10);
        // normal should be perpendicular to both dPdu and dPdv
        assert!(dpu.dot(p.normal_at(0.0, 0.0)).abs() < 1e-10);
        assert!(dpv.dot(p.normal_at(0.0, 0.0)).abs() < 1e-10);
        // point's Z should be 0
        assert!(pt.z.abs() < 1e-10);
    }

    // ── Cylinder ────────────────────────────────────────────────────────

    #[test]
    fn cylinder_derivatives_and_normal() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        let (p, dpu, dpv) = cyl.derivatives(0.0, 5.0);
        // dP/dv should be axis (Z)
        assert!((dpv - DVec3::Z).length() < 1e-10);
        // dP/du should be tangent around the cylinder, perpendicular to radius
        let radial = p - DVec3::new(0.0, 0.0, 5.0);
        let radial = radial.normalize_or_zero();
        assert!(dpu.dot(radial).abs() < 1e-10);
        // normal should be the radial direction
        let n = cyl.normal_at(0.0, 5.0);
        assert!((n - radial).length() < 1e-10);
    }

    #[test]
    fn cylinder_normal_radial() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        // Cylinder with axis Z, ref_dir=X:
        // x_ax = X, y_ax = Z×X = Y
        // At u=0, normal should point in the X direction (radial outward)
        let n = cyl.normal_at(0.0, 0.0);
        assert!((n - DVec3::X).length() < 1e-10);
    }

    // ── Sphere ──────────────────────────────────────────────────────────

    #[test]
    fn sphere_derivatives_and_normal() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        // At equator (v=PI/2), u=0: point = (5, 0, 0)
        use std::f64::consts::PI;
        let (p, dpu, dpv) = s.derivatives(0.0, PI / 2.0);
        // Point radius should be 5
        assert!((p.length() - 5.0).abs() < 1e-10);
        // dP/du should be perpendicular to point
        assert!(dpu.dot(p.normalize_or_zero()).abs() < 1e-10);
        // dP/dv at equator should point toward -Z
        assert!((dpv - DVec3::new(0.0, 0.0, -5.0)).length() < 1e-10);
        // Normal should point radially outward
        let n = s.normal_at(0.0, PI / 2.0);
        assert!((n - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn sphere_normal_at_poles() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        // North pole v=0: point=(0,0,5), normal should be +Z
        let n_north = s.normal_at(0.0, 0.0);
        assert!((n_north - DVec3::Z).length() < 1e-10);
        // South pole v=PI: point=(0,0,-5), normal should be -Z
        let n_south = s.normal_at(0.0, std::f64::consts::PI);
        assert!((n_south + DVec3::Z).length() < 1e-10);
    }

    // ── Cone ────────────────────────────────────────────────────────────

    #[test]
    fn cone_derivatives() {
        let sa = 30.0_f64.to_radians(); // 30 degree half-angle
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: sa,
        };
        let (_p, dpu, dpv) = cone.derivatives(0.0, 0.0);
        // At v=0: radial = 2, axial = 0
        // dP/du at u=0: radial * (-sin(0)*x_ax + cos(0)*y_ax)
        //   where x_ax = any_perpendicular(Z) = Y, y_ax = Z×Y = -X
        //   = 2 * (0*Y + 1*(-X)) = 2*(-X)
        // dP/dv at v=0: da*axis + dr*r_vec = cos(sa)*Z + sin(sa)*Y
        //   where r_vec = cos(0)*Y + sin(0)*(-X) = Y
        //   = cos(sa)*Z + sin(sa)*Y
        assert!(dpu.length() > 0.0);
        assert!(dpv.length() > 0.0);
        // dP/du should be perpendicular to the radial direction
        let n = cone.normal_at(0.0, 0.0);
        assert!(dpu.dot(n).abs() < 1e-10);
        assert!(dpv.dot(n).abs() < 1e-10);
    }

    // ── Torus ───────────────────────────────────────────────────────────

    #[test]
    fn torus_derivatives_and_normal() {
        use std::f64::consts::PI;
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // OCCT-aligned: x_ax=X, y_ax=Y. At u=0, v=0: outer equator, point=(6,0,0)
        let (p, dpu, dpv) = t.derivatives(0.0, 0.0);
        assert!((p - DVec3::new(6.0, 0.0, 0.0)).length() < 1e-9);
        // dP/du should be in the Y direction (major circle tangent)
        assert!((dpu - DVec3::new(0.0, 6.0, 0.0)).length() < 1e-9);
        // dP/dv should be in the Z direction (minor circle tangent at v=0)
        assert!((dpv - DVec3::Z).length() < 1e-9);
        // Normal at outer equator should be outward radial (X)
        let n = t.normal_at(0.0, 0.0);
        assert!((n - DVec3::X).length() < 1e-9);
    }

    #[test]
    fn torus_inner_equator_normal() {
        use std::f64::consts::PI;
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Inner equator: u=0, v=PI → point=(4,0,0), normal should be -X (inward)
        let n_inner = t.normal_at(0.0, PI);
        assert!((n_inner + DVec3::X).length() < 1e-9);
    }

    // ── BSplineSurface ──────────────────────────────────────────────────

    #[test]
    fn bspline_surface_eval_d0() {
        let surf = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 10.0, 0.0)],
                vec![DVec3::new(10.0, 0.0, 0.0), DVec3::new(10.0, 10.0, 0.0)],
            ],
            weights: vec![vec![1.0; 2]; 2],
        };
        // Degree 1 surface = bilinear, should interpolate corners
        assert!((surf.point_at(0.0, 0.0) - DVec3::ZERO).length() < 1e-10);
        assert!((surf.point_at(1.0, 1.0) - DVec3::new(10.0, 10.0, 0.0)).length() < 1e-10);
        assert!((surf.point_at(0.5, 0.5) - DVec3::new(5.0, 5.0, 0.0)).length() < 1e-10);
    }

    // ── Surface3 dispatch verification ───────────────────────────────────

    #[test]
    fn surface3_plane_dispatch() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        // OCCT-aligned: normal=Z gives u_dir=X, v_dir=Y
        let (_p, dpu, dpv) = s.derivatives(1.0, 2.0);
        assert!((dpu - DVec3::X).length() < 1e-10);
        assert!((dpv - DVec3::Y).length() < 1e-10);
        assert!(s.normal_at(1.0, 2.0) == DVec3::Z);
    }

    #[test]
    fn surface3_sphere_dispatch() {
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        });
        // Verify D1/dPdu through Surface3 dispatch
        let (p, _dpu, _dpv) = s.derivatives(0.0, 0.0);
        assert!((p - DVec3::new(0.0, 0.0, 3.0)).length() < 1e-9);
        let n = s.normal_at(0.0, 0.0);
        assert!((n - DVec3::Z).length() < 1e-9);
    }

    #[test]
    fn surface3_cylinder_dispatch() {
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: DVec3::X,
        });
        // Cylinder with axis Z, ref_dir=X: dP/du at u=0 = R*y_ax = 2*Y = (0,2,0)
        let (_p, dpu, dpv) = s.derivatives(0.0, 0.0);
        // dP/dv should be axis direction
        assert!((dpv - DVec3::Z).length() < 1e-10);
        // dP/du should be tangent (perpendicular to radial)
        assert!((dpu - DVec3::new(0.0, 2.0, 0.0)).length() < 1e-10);
    }
}
