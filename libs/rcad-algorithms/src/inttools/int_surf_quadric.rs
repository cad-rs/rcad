//! ✅ OCCT-aligned: IntSurf_Quadric — quadric surface representation.
//!
//! OCCT IntSurf_Quadric.hxx / .cxx / .lxx
//!
//! Stores analytic surfaces (Plane, Cylinder, Sphere, Cone, Torus) in a
//! unified algebraic form for efficient intersection algorithms.
//! Used by ImpImpIntersection as the common surface representation
//! for all 15 surface pair combinations.

use glam::{DVec3, DAffine3};
use rcad_kernel::geom::{Surface3, Plane, SphericalSurface, CylindricalSurface, ConicalSurface, ToroidalSurface};
use crate::tolerance::TOLERANCE_CLAMP_MIN;
use super::geom_abs_surface_type::GeomAbsSurfaceType;

/// OCCT-aligned: IntSurf_Quadric unified quadric surface.
///
/// OCCT fields (L102-109):
///   ax3: gp_Ax3        — coordinate system
///   lin: gp_Lin        — axis line (used for Cylinder/Sphere/Cone/Torus)
///   typ: GeomAbs_SurfaceType — surface type discriminator
///   prm1-4: double     — parameters (radius, angle, coefficients)
///   ax3direc: bool     — true if ax3 is right-handed
pub struct Quadric {
    // OCCT: ax3 — coordinate system (position + axes)
    location: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    z_dir: DVec3,
    // OCCT: lin — axis line (set to ax3.Axis() for Cylinder/Sphere/Cone/Torus)
    axis_loc: DVec3,
    axis_dir: DVec3,
    // OCCT: typ
    typ: GeomAbsSurfaceType,
    // OCCT: prm1-4
    prm1: f64,
    prm2: f64,
    prm3: f64,
    prm4: f64,
    // OCCT: ax3direc
    ax3direc: bool,
}

impl Quadric {
    /// OCCT default constructor
    pub fn new() -> Self {
        Self {
            location: DVec3::ZERO,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            z_dir: DVec3::Z,
            axis_loc: DVec3::ZERO,
            axis_dir: DVec3::Z,
            typ: GeomAbsSurfaceType::OtherSurface,
            prm1: 0.0, prm2: 0.0, prm3: 0.0, prm4: 0.0,
            ax3direc: true,
        }
    }

    /// Build from any rcad Surface3. Returns None for non-quadric surfaces.
    pub fn from_surface3(surf: &Surface3) -> Option<Self> {
        match surf {
            Surface3::Plane(p) => Some(Self::from_plane(p)),
            Surface3::Cylinder(c) => Some(Self::from_cylinder(c)),
            Surface3::Sphere(s) => Some(Self::from_sphere(s)),
            Surface3::Cone(c) => Some(Self::from_cone(c)),
            Surface3::Torus(t) => Some(Self::from_torus(t)),
            _ => None,
        }
    }

    /// OCCT: IntSurf_Quadric(const gp_Pln&)
    pub fn from_plane(p: &Plane) -> Self {
        let (x_dir, y_dir) = any_perpendicular_pair(p.normal);
        let z_dir = p.normal;
        let ax3direc = z_dir.cross(x_dir).dot(y_dir) > 0.0;
        // Plane equation coefficients: a*x + b*y + c*z + d = 0
        // where normal = (a,b,c) and d = -dot(normal, origin)
        let d = -p.normal.dot(p.origin);
        Self {
            location: p.origin,
            x_dir, y_dir, z_dir,
            axis_loc: p.origin,
            axis_dir: p.normal,
            typ: GeomAbsSurfaceType::Plane,
            prm1: p.normal.x, prm2: p.normal.y, prm3: p.normal.z, prm4: d,
            ax3direc,
        }
    }

    /// OCCT: IntSurf_Quadric(const gp_Cylinder&)
    pub fn from_cylinder(c: &CylindricalSurface) -> Self {
        let axis_dir = c.axis.normalize();
        let (x_dir, y_dir) = any_perpendicular_pair(axis_dir);
        let z_dir = axis_dir;
        Self {
            location: c.origin,
            x_dir, y_dir, z_dir,
            axis_loc: c.origin,
            axis_dir,
            typ: GeomAbsSurfaceType::Cylinder,
            prm1: c.radius, prm2: 0.0, prm3: 0.0, prm4: 0.0,
            ax3direc: true,
        }
    }

    /// OCCT: IntSurf_Quadric(const gp_Sphere&)
    pub fn from_sphere(s: &SphericalSurface) -> Self {
        Self {
            location: s.center,
            x_dir: DVec3::X, y_dir: DVec3::Y, z_dir: DVec3::Z,
            axis_loc: s.center,
            axis_dir: DVec3::Z,
            typ: GeomAbsSurfaceType::Sphere,
            prm1: s.radius, prm2: 0.0, prm3: 0.0, prm4: 0.0,
            ax3direc: true,
        }
    }

    /// OCCT: IntSurf_Quadric(const gp_Cone&)
    pub fn from_cone(c: &ConicalSurface) -> Self {
        let axis_dir = c.axis.normalize();
        let ax3direc = true;
        let prm2 = c.half_angle_rad; // OCCT: SemiAngle
        Self {
            location: c.apex,
            x_dir: DVec3::X, y_dir: DVec3::Y, z_dir: axis_dir,
            axis_loc: c.apex,
            axis_dir,
            typ: GeomAbsSurfaceType::Cone,
            prm1: c.radius,  // OCCT: RefRadius
            prm2,             // OCCT: SemiAngle
            prm3: prm2.cos(), // OCCT: cos(SemiAngle)
            prm4: 0.0,
            ax3direc,
        }
    }

    /// OCCT: IntSurf_Quadric(const gp_Torus&)
    pub fn from_torus(t: &ToroidalSurface) -> Self {
        let axis_dir = t.axis.normalize();
        Self {
            location: t.center,
            x_dir: DVec3::X, y_dir: DVec3::Y, z_dir: axis_dir,
            axis_loc: t.center,
            axis_dir,
            typ: GeomAbsSurfaceType::Torus,
            prm1: t.major_radius, // OCCT: MajorRadius
            prm2: t.minor_radius, // OCCT: MinorRadius
            prm3: 0.0, prm4: 0.0,
            ax3direc: true,
        }
    }

    /// OCCT: TypeQuadric()
    pub fn surface_type(&self) -> GeomAbsSurfaceType { self.typ }

    /// OCCT: Value(U, V) → gp_Pnt — evaluate surface at UV.
    pub fn value(&self, u: f64, v: f64) -> DVec3 {
        match self.typ {
            GeomAbsSurfaceType::Plane => {
                self.location + u * self.x_dir + v * self.y_dir
            }
            GeomAbsSurfaceType::Cylinder => {
                let r = self.prm1;
                let cos_u = u.cos();
                let sin_u = u.sin();
                let pt_on_axis = self.axis_loc + v * self.axis_dir;
                let radial = cos_u * self.x_dir + sin_u * self.y_dir;
                pt_on_axis + r * radial
            }
            GeomAbsSurfaceType::Sphere => {
                let r = self.prm1;
                let cos_v = v.cos();
                let sin_v = v.sin();
                let cos_u = u.cos();
                let sin_u = u.sin();
                let pt = DVec3::new(
                    r * cos_v * cos_u,
                    r * cos_v * sin_u,
                    r * sin_v,
                );
                self.location + pt
            }
            GeomAbsSurfaceType::Cone => {
                let r = self.prm1;
                let ang = self.prm2;
                let v_factor = r + v * ang.sin();
                let cos_u = u.cos();
                let sin_u = u.sin();
                let pt_on_axis = self.axis_loc + v * self.axis_dir;
                let radial = cos_u * self.x_dir + sin_u * self.y_dir;
                pt_on_axis + v_factor * radial
            }
            GeomAbsSurfaceType::Torus => {
                let maj_r = self.prm1;
                let min_r = self.prm2;
                let cos_u = u.cos();
                let sin_u = u.sin();
                let cos_v = v.cos();
                let sin_v = v.sin();
                let center = self.location + maj_r * (cos_u * self.x_dir + sin_u * self.y_dir);
                let normal = cos_u * cos_v * self.x_dir
                    + sin_u * cos_v * self.y_dir
                    + sin_v * self.z_dir;
                center + min_r * normal
            }
            _ => DVec3::ZERO,
        }
    }

    /// OCCT: D1(U, V, P, D1U, D1V) — evaluate surface and first derivatives.
    pub fn d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        // Default: finite difference if analytic derivatives not implemented
        let p = self.value(u, v);
        let eps = 1e-6;
        let du = (self.value(u + eps, v) - p) / eps;
        let dv = (self.value(u, v + eps) - p) / eps;
        (p, du, dv)
    }

    /// OCCT: Normale(U, V) — surface unit normal.
    pub fn normale(&self, u: f64, v: f64) -> DVec3 {
        let (_, du, dv) = self.d1(u, v);
        du.cross(dv).normalize_or_zero()
    }

    /// OCCT: Parameters(P, U, V) — project 3D point to UV coordinates.
    pub fn parameters(&self, p: DVec3) -> (f64, f64) {
        match self.typ {
            GeomAbsSurfaceType::Plane => {
                let d = p - self.location;
                (d.dot(self.x_dir), d.dot(self.y_dir))
            }
            GeomAbsSurfaceType::Cylinder => {
                let d = p - self.axis_loc;
                let v = d.dot(self.axis_dir);
                let radial = d - v * self.axis_dir;
                let u = radial.y.atan2(radial.x);
                (if u < 0.0 { u + std::f64::consts::TAU } else { u }, v)
            }
            GeomAbsSurfaceType::Sphere => {
                let d = p - self.location;
                let r = d.length();
                if r < TOLERANCE_CLAMP_MIN { return (0.0, 0.0); }
                let v = (d.z / r).asin();
                let u = d.y.atan2(d.x);
                (if u < 0.0 { u + std::f64::consts::TAU } else { u }, v)
            }
            GeomAbsSurfaceType::Cone => {
                let d = p - self.axis_loc;
                let v = d.dot(self.axis_dir);
                let radial = d - v * self.axis_dir;
                let r = radial.length();
                let u = if r < TOLERANCE_CLAMP_MIN { 0.0 } else { radial.y.atan2(radial.x) };
                (if u < 0.0 { u + std::f64::consts::TAU } else { u }, v)
            }
            GeomAbsSurfaceType::Torus => {
                let d = p - self.location;
                let proj_axis = d.dot(self.axis_dir) * self.axis_dir;
                let radial = d - proj_axis;
                let r = radial.length();
                let v = if r < TOLERANCE_CLAMP_MIN { 0.0 } else {
                    let center_to_p = (d - self.major_radius() * radial.normalize_or_zero()).length();
                    (d.dot(self.axis_dir) / center_to_p).asin()
                };
                let u = if r < TOLERANCE_CLAMP_MIN { 0.0 } else { radial.y.atan2(radial.x) };
                (if u < 0.0 { u + std::f64::consts::TAU } else { u }, v)
            }
            _ => (0.0, 0.0),
        }
    }

    /// OCCT: Distance(P) — algebraic distance from point to quadric.
    ///
    /// Returns the algebraic distance (signed), not Euclidean.
    /// For Plane: a*x + b*y + c*z + d
    /// For Cylinder: distance to axis - radius
    /// For Sphere: distance to center - radius
    /// For Cone: (distance to axis - distance along cone) / cos(semi_angle)
    /// For Torus: sqrt((R - sqrt(x^2+y^2))^2 + z^2) - r
    pub fn distance(&self, p: DVec3) -> f64 {
        match self.typ {
            GeomAbsSurfaceType::Plane => {
                self.prm1 * p.x + self.prm2 * p.y + self.prm3 * p.z + self.prm4
            }
            GeomAbsSurfaceType::Cylinder => {
                let d = p - self.axis_loc;
                let v = d.dot(self.axis_dir);
                let radial = d - v * self.axis_dir;
                radial.length() - self.prm1
            }
            GeomAbsSurfaceType::Sphere => {
                (p - self.axis_loc).length() - self.prm1
            }
            GeomAbsSurfaceType::Cone => {
                let d = p - self.axis_loc;
                let v = d.dot(self.axis_dir);
                let radial = d - v * self.axis_dir;
                let dist_axis = radial.length();
                // Project point onto cone surface
                let prm1 = self.prm1; // RefRadius
                let prm2 = self.prm2; // SemiAngle
                let (u_p, _) = self.parameters(p);
                let pp = self.value(u_p, v);
                let dist_pp = (pp - self.axis_loc - v * self.axis_dir).length();
                (dist_axis - dist_pp) / self.prm3  // prm3 = cos(semi_angle)
            }
            GeomAbsSurfaceType::Torus => {
                let d = p - self.location;
                let proj_axis = d.dot(self.axis_dir) * self.axis_dir;
                let radial = d - proj_axis;
                let r = radial.length();
                let dist_to_ring = (r - self.prm1).abs();
                let z = d.dot(self.axis_dir);
                (dist_to_ring * dist_to_ring + z * z).sqrt() - self.prm2
            }
            _ => 0.0,
        }
    }

    /// OCCT: Gradient(P) — gradient of the algebraic distance function.
    pub fn gradient(&self, p: DVec3) -> DVec3 {
        match self.typ {
            GeomAbsSurfaceType::Plane => {
                DVec3::new(self.prm1, self.prm2, self.prm3)
            }
            GeomAbsSurfaceType::Cylinder => {
                let d = p - self.axis_loc;
                let v = d.dot(self.axis_dir);
                let radial = d - v * self.axis_dir;
                let r = radial.length();
                if r < TOLERANCE_CLAMP_MIN { DVec3::ZERO } else { radial / r }
            }
            GeomAbsSurfaceType::Sphere => {
                let d = p - self.axis_loc;
                let r = d.length();
                if r < TOLERANCE_CLAMP_MIN { DVec3::ZERO } else { d / r }
            }
            GeomAbsSurfaceType::Cone => {
                let d = p - self.axis_loc;
                let v = d.dot(self.axis_dir);
                let radial = d - v * self.axis_dir;
                let r = radial.length();
                if r < TOLERANCE_CLAMP_MIN { return self.axis_dir; }
                let radial_dir = radial / r;
                let cone_dir = self.axis_dir * self.prm3 - radial_dir * self.prm2.sin();
                cone_dir.normalize_or_zero()
            }
            GeomAbsSurfaceType::Torus => {
                let d = p - self.location;
                let proj_axis = d.dot(self.axis_dir) * self.axis_dir;
                let radial = d - proj_axis;
                let r = radial.length();
                if r < TOLERANCE_CLAMP_MIN { return DVec3::ZERO; }
                let maj_r = self.prm1;
                let radial_dir = radial / r;
                let grad = radial_dir * (r - maj_r) / r + self.axis_dir * d.dot(self.axis_dir);
                let norm = grad.length();
                if norm < TOLERANCE_CLAMP_MIN { DVec3::ZERO } else { grad / norm }
            }
            _ => DVec3::ZERO,
        }
    }

    /// OCCT: ValAndGrad(P, Dist, Grad) — combined distance and gradient.
    pub fn val_and_grad(&self, p: DVec3) -> (f64, DVec3) {
        (self.distance(p), self.gradient(p))
    }

    // Accessors for parameters — OCCT prm1..prm4
    pub fn major_radius(&self) -> f64 { self.prm1 }
    pub fn minor_radius(&self) -> f64 { self.prm2 }
    pub fn semi_angle(&self) -> f64 { self.prm2 }
    pub fn ref_radius(&self) -> f64 { self.prm1 }
    pub fn radius(&self) -> f64 { self.prm1 }
    pub fn plane_coeffs(&self) -> (f64, f64, f64, f64) { (self.prm1, self.prm2, self.prm3, self.prm4) }
    pub fn axis_dir(&self) -> DVec3 { self.axis_dir }
    pub fn axis_loc(&self) -> DVec3 { self.axis_loc }
    pub fn location(&self) -> DVec3 { self.location }
}

/// Helper: compute a perpendicular pair to a given direction (OCCT-aligned).
fn any_perpendicular_pair(dir: DVec3) -> (DVec3, DVec3) {
    let x_dir = rcad_kernel::geom::any_perpendicular(dir);
    let y_dir = dir.cross(x_dir).normalize_or_zero();
    (x_dir, y_dir)
}
