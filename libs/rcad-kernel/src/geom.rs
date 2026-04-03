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
    pub direction: Vec3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Circle3 {
    pub center: Point3,
    pub normal: Vec3,
    pub radius: f64,
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Curve3 {
    Line(Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
    BSpline(BSplineCurve3),
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Plane {
    pub origin: Point3,
    pub normal: Vec3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CylindricalSurface {
    pub origin: Point3,
    pub axis: Vec3,
    pub radius: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SphericalSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ConicalSurface {
    pub apex: Point3,
    pub axis: Vec3,
    pub radius: f64,
    pub half_angle_rad: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ToroidalSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub major_radius: f64,
    pub minor_radius: f64,
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Surface3 {
    Plane(Plane),
    Cylinder(CylindricalSurface),
    Sphere(SphericalSurface),
    Cone(ConicalSurface),
    Torus(ToroidalSurface),
    BSpline(BSplineSurface),
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

// ── 2D Geometry (parameter-space / PCurve types) ─────────────────────────────

/// A line in 2D parameter space: point + direction.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Line2d {
    pub origin: Point2,
    pub direction: Vec2,
}

/// A circle in 2D parameter space.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Circle2d {
    pub center: Point2,
    pub radius: f64,
}

/// A curve defined in the 2D parameter space (u, v) of a surface.
///
/// Used for PCurves: the image of a 3D edge on the parameter domain of an
/// adjacent face surface. Analogous to OCCT `Geom2d_Curve`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Curve2d {
    Line(Line2d),
    Circle(Circle2d),
}

// ── Geometric evaluation traits ──────────────────────────────────────────────

/// Returns a vector perpendicular to `v`. Stable for any non-zero input.
pub fn any_perpendicular(v: DVec3) -> DVec3 {
    // Pick the axis least aligned with v, then cross.
    let abs = v.abs();
    let candidate = if abs.x <= abs.y && abs.x <= abs.z {
        DVec3::X
    } else if abs.y <= abs.z {
        DVec3::Y
    } else {
        DVec3::Z
    };
    v.cross(candidate).normalize()
}

/// Parametric evaluation of a 3D curve: `t → Point3`.
///
/// Mirrors OCCT `Geom_Curve::Value(t)` / `D1(t)`.
pub trait CurveEval {
    /// Point on the curve at parameter `t`.
    fn point_at(&self, t: f64) -> DVec3;
    /// Unit tangent vector at parameter `t`.
    fn tangent_at(&self, t: f64) -> DVec3;
    /// Natural parameter domain `[t_min, t_max]`.
    /// Lines use `[NEG_INFINITY, INFINITY]`; circles/ellipses use `[0, 2π]`.
    fn default_domain(&self) -> [f64; 2];
}

/// Parametric evaluation of a 3D surface: `(u, v) → Point3`.
///
/// Mirrors OCCT `Geom_Surface::Value(u, v)`.
pub trait SurfaceEval {
    /// Point on the surface at parameter `(u, v)`.
    fn point_at(&self, u: f64, v: f64) -> DVec3;
    /// Outward unit normal at parameter `(u, v)`.
    fn normal_at(&self, u: f64, v: f64) -> DVec3;
    /// Natural parameter domain `[u_min, u_max, v_min, v_max]`.
    fn default_domain(&self) -> [f64; 4];
}

/// Parametric evaluation of a 2D curve (PCurve): `t → Point2`.
pub trait Curve2dEval {
    fn point_at(&self, t: f64) -> DVec2;
}

// ── CurveEval implementations ─────────────────────────────────────────────────

impl CurveEval for Line3 {
    fn point_at(&self, t: f64) -> DVec3 {
        self.origin + t * self.direction
    }
    fn tangent_at(&self, _t: f64) -> DVec3 {
        self.direction
    }
    fn default_domain(&self) -> [f64; 2] {
        [f64::NEG_INFINITY, f64::INFINITY]
    }
}

impl CurveEval for Circle3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.normal);
        let y_ax = self.normal.cross(x_ax);
        self.center + self.radius * (t.cos() * x_ax + t.sin() * y_ax)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.normal);
        let y_ax = self.normal.cross(x_ax);
        (-t.sin() * x_ax + t.cos() * y_ax).normalize()
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
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
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
    }
}

impl CurveEval for Curve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.point_at(t),
            Curve3::Circle(c) => c.point_at(t),
            Curve3::Ellipse(c) => c.point_at(t),
            Curve3::BSpline(c) => c.point_at(t),
        }
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.tangent_at(t),
            Curve3::Circle(c) => c.tangent_at(t),
            Curve3::Ellipse(c) => c.tangent_at(t),
            Curve3::BSpline(c) => c.tangent_at(t),
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        match self {
            Curve3::Line(c) => c.default_domain(),
            Curve3::Circle(c) => c.default_domain(),
            Curve3::Ellipse(c) => c.default_domain(),
            Curve3::BSpline(c) => c.default_domain(),
        }
    }
}

// ── SurfaceEval implementations ───────────────────────────────────────────────

impl SurfaceEval for Plane {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.normal);
        let y_ax = self.normal.cross(x_ax);
        self.origin + u * x_ax + v * y_ax
    }
    fn normal_at(&self, _u: f64, _v: f64) -> DVec3 {
        self.normal
    }
    fn default_domain(&self) -> [f64; 4] {
        [f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY]
    }
}

impl SurfaceEval for CylindricalSurface {
    /// u = azimuth angle [0, 2π], v = height along axis.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        self.origin + self.radius * (u.cos() * x_ax + u.sin() * y_ax) + v * self.axis
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        (u.cos() * x_ax + u.sin() * y_ax).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, f64::NEG_INFINITY, f64::INFINITY]
    }
}

impl SurfaceEval for SphericalSurface {
    /// u = longitude [0, 2π], v = colatitude [0, π] (0 = north pole).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        self.center
            + self.radius
                * (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * self.axis)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * self.axis).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, PI]
    }
}

impl SurfaceEval for ConicalSurface {
    /// u = azimuth [0, 2π], v = distance along slant from apex (v ≥ 0).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        let r = v * self.half_angle_rad.tan();
        self.apex + v * self.axis + r * (u.cos() * x_ax + u.sin() * y_ax)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        let radial = u.cos() * x_ax + u.sin() * y_ax;
        let half = self.half_angle_rad;
        (radial * half.cos() - self.axis * half.sin()).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, f64::INFINITY]
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
}

impl SurfaceEval for Surface3 {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface3::Plane(s) => s.point_at(u, v),
            Surface3::Cylinder(s) => s.point_at(u, v),
            Surface3::Sphere(s) => s.point_at(u, v),
            Surface3::Cone(s) => s.point_at(u, v),
            Surface3::Torus(s) => s.point_at(u, v),
            Surface3::BSpline(s) => s.point_at(u, v),
        }
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface3::Plane(s) => s.normal_at(u, v),
            Surface3::Cylinder(s) => s.normal_at(u, v),
            Surface3::Sphere(s) => s.normal_at(u, v),
            Surface3::Cone(s) => s.normal_at(u, v),
            Surface3::Torus(s) => s.normal_at(u, v),
            Surface3::BSpline(s) => s.normal_at(u, v),
        }
    }
    fn default_domain(&self) -> [f64; 4] {
        match self {
            Surface3::Plane(s) => s.default_domain(),
            Surface3::Cylinder(s) => s.default_domain(),
            Surface3::Sphere(s) => s.default_domain(),
            Surface3::Cone(s) => s.default_domain(),
            Surface3::Torus(s) => s.default_domain(),
            Surface3::BSpline(s) => s.default_domain(),
        }
    }
}

// ── BSpline evaluation ────────────────────────────────────────────────────────

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
        for i in degree..knots.len() - degree - 1 {
            if knots[i] <= t_clamped {
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
            for c in 0..4 {
                d[j][c] = (1.0 - alpha) * d[j - 1][c] + alpha * d[j][c];
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

impl CurveEval for BSplineCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        de_boor(self.degree, &self.knots, &self.control_points, &self.weights, t)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        // Central difference approximation
        let eps = (self.default_domain()[1] - self.default_domain()[0]) * 1e-6;
        let t0 = self.default_domain()[0];
        let t1 = self.default_domain()[1];
        let t_lo = (t - eps).max(t0);
        let t_hi = (t + eps).min(t1);
        let dp = self.point_at(t_hi) - self.point_at(t_lo);
        let len = dp.length();
        if len < 1e-15 { DVec3::X } else { dp / len }
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
        // Tensor product: evaluate in u for each v row, then interpolate in v
        let n_u = self.control_points.len();
        if n_u == 0 {
            return DVec3::ZERO;
        }
        let n_v = self.control_points[0].len();
        if n_v == 0 {
            return DVec3::ZERO;
        }
        // Evaluate each v-row of control points along u
        let row_points: Vec<DVec3> = (0..n_v)
            .map(|j| {
                let pts: Vec<DVec3> = (0..n_u).map(|i| self.control_points[i][j]).collect();
                let wts: Vec<f64> = (0..n_u).map(|i| self.weights[i][j]).collect();
                de_boor(self.degree_u, &self.knots_u, &pts, &wts, u)
            })
            .collect();
        let unit_weights = vec![1.0; n_v];
        de_boor(self.degree_v, &self.knots_v, &row_points, &unit_weights, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let [u0, u1, v0, v1] = self.default_domain();
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
        let u1 = if nu > du + 1 { self.knots_u[nu - du - 1] } else { 1.0 };
        let v0 = if nv > dv { self.knots_v[dv] } else { 0.0 };
        let v1 = if nv > dv + 1 { self.knots_v[nv - dv - 1] } else { 1.0 };
        [u0, u1, v0, v1]
    }
}

// ── Curve2dEval implementations ───────────────────────────────────────────────

impl Curve2dEval for Line2d {
    fn point_at(&self, t: f64) -> DVec2 {
        self.origin + t * self.direction
    }
}

impl Curve2dEval for Circle2d {
    fn point_at(&self, t: f64) -> DVec2 {
        self.center + self.radius * DVec2::new(t.cos(), t.sin())
    }
}

impl Curve2dEval for Curve2d {
    fn point_at(&self, t: f64) -> DVec2 {
        match self {
            Curve2d::Line(c) => c.point_at(t),
            Curve2d::Circle(c) => c.point_at(t),
        }
    }
}

#[cfg(test)]
mod eval_tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn line3_point_at() {
        let l = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
        assert!((l.point_at(3.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn circle3_point_at_zero_is_on_circle() {
        // Circle in XY plane, normal = Z
        let c = Circle3 { center: DVec3::ZERO, normal: DVec3::Z, radius: 2.0 };
        let p0 = c.point_at(0.0);
        assert!((p0.length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn circle3_full_revolution_closes() {
        let c = Circle3 { center: DVec3::new(1.0, 2.0, 3.0), normal: DVec3::Y, radius: 5.0 };
        let p0 = c.point_at(0.0);
        let p2pi = c.point_at(2.0 * PI);
        assert!((p0 - p2pi).length() < 1e-10);
    }

    #[test]
    fn circle3_quarter_turn() {
        let c = Circle3 { center: DVec3::ZERO, normal: DVec3::Z, radius: 1.0 };
        let p0 = c.point_at(0.0);
        let p90 = c.point_at(FRAC_PI_2);
        // 90° rotation: p0 and p90 should be perpendicular from center
        assert!((p0.dot(p90)).abs() < 1e-10);
        assert!((p90.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_surface_north_pole() {
        let s = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Y, radius: 3.0 };
        // v=0 is north pole regardless of u
        let p = s.point_at(0.0, 0.0);
        // Should be at (0, 3, 0)
        assert!((p - DVec3::new(0.0, 3.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn sphere_surface_point_on_sphere() {
        let s = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Y, radius: 2.0 };
        for u in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            for v in [0.1, 0.5, 1.0, PI / 2.0, PI - 0.1] {
                let p = s.point_at(u, v);
                assert!((p.length() - 2.0).abs() < 1e-9, "u={u} v={v} |p|={}", p.length());
            }
        }
    }

    #[test]
    fn cylinder_surface_point_on_cylinder() {
        let c = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Y, radius: 3.0 };
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
        let t = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Y, major_radius: 5.0, minor_radius: 1.0 };
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
}
