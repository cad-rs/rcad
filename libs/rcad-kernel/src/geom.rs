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

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Curve3 {
    Line(Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Surface3 {
    Plane(Plane),
    Cylinder(CylindricalSurface),
    Sphere(SphericalSurface),
    Cone(ConicalSurface),
    Torus(ToroidalSurface),
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
        }
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.tangent_at(t),
            Curve3::Circle(c) => c.tangent_at(t),
            Curve3::Ellipse(c) => c.tangent_at(t),
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        match self {
            Curve3::Line(c) => c.default_domain(),
            Curve3::Circle(c) => c.default_domain(),
            Curve3::Ellipse(c) => c.default_domain(),
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
        }
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface3::Plane(s) => s.normal_at(u, v),
            Surface3::Cylinder(s) => s.normal_at(u, v),
            Surface3::Sphere(s) => s.normal_at(u, v),
            Surface3::Cone(s) => s.normal_at(u, v),
            Surface3::Torus(s) => s.normal_at(u, v),
        }
    }
    fn default_domain(&self) -> [f64; 4] {
        match self {
            Surface3::Plane(s) => s.default_domain(),
            Surface3::Cylinder(s) => s.default_domain(),
            Surface3::Sphere(s) => s.default_domain(),
            Surface3::Cone(s) => s.default_domain(),
            Surface3::Torus(s) => s.default_domain(),
        }
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
