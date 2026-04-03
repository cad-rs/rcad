use glam::DVec3;
use serde::{Deserialize, Serialize};

pub type Point3 = DVec3;
pub type Vec3 = DVec3;

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
        u_segments: usize,
        v_segments: usize,
    },
    Cylinder {
        radius: f64,
        height: f64,
        segments: usize,
    },
    Cone {
        base_radius: f64,
        height: f64,
        segments: usize,
    },
    Torus {
        major_radius: f64,
        minor_radius: f64,
        major_segments: usize,
        minor_segments: usize,
    },
}
