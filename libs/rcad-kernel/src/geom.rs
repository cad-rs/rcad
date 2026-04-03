use glam::DVec3;
use serde::{Deserialize, Serialize};

pub type Point3 = DVec3;
pub type Vec3 = DVec3;
pub type Point2 = glam::DVec2;
pub type Vec2 = glam::DVec2;

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
