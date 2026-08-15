//! 3D curve projection onto analytic surfaces (ProjLib).
//!
//! Projects elementary curves (line, circle, ellipse, parabola, hyperbola) onto
//! analytic surfaces (plane, cylinder, sphere, cone, torus).
//!
//! OCCT TKGeomBase ProjLib package: ProjLib_Projector, ProjLib_Plane,
//! ProjLib_Cylinder, ProjLib_Sphere, ProjLib_Cone, ProjLib_Torus.

#![allow(clippy::manual_clamp)]

use glam::{DVec2, DVec3};

use crate::geom::{
    Circle3, ConicalSurface, Curve2d, Curve3, CylindricalSurface, Ellipse3, Hyperbola3, Line3,
    Parabola3, Plane, SphericalSurface, Surface3, ToroidalSurface,
};

// ============================================================================
// CurveType — mirrors GeomAbs_CurveType
// ============================================================================

/// Type of the resulting projected curve.
///
/// OCCT: `GeomAbs_CurveType`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurveType {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    Bezier,
    BSpline,
    Other,
}

// ============================================================================
// Projector — base class storing projection result
// ============================================================================

/// Base class for projection algorithms.
///
/// OCCT: `ProjLib_Projector`.
pub struct Projector {
    proj_type: CurveType,
    lin: Line3,
    circ: Circle3,
    elips: Ellipse3,
    hypr: Hyperbola3,
    parab: Parabola3,
    bspline: Option<Curve2d>,
    is_periodic: bool,
    is_done: bool,
}

impl Projector {
    /// Create an empty projector.
    ///
    /// OCCT: default constructor.
    pub fn new() -> Self {
        Projector {
            proj_type: CurveType::Other,
            lin: Line3::new(DVec3::ZERO, DVec3::X),
            circ: Circle3::new(DVec3::ZERO, DVec3::Z, 1.0),
            elips: Ellipse3 {
                center: DVec3::ZERO,
                normal: DVec3::Z,
                major_dir: DVec3::X,
                major_radius: 2.0,
                minor_radius: 1.0,
            },
            hypr: Hyperbola3 {
                center: DVec3::ZERO,
                normal: DVec3::Z,
                major_dir: DVec3::X,
                semi_major: 2.0,
                semi_minor: 1.0,
            },
            parab: Parabola3 {
                vertex: DVec3::ZERO,
                normal: DVec3::Z,
                axis_dir: DVec3::X,
                focal_param: 1.0,
            },
            bspline: None,
            is_periodic: false,
            is_done: false,
        }
    }

    /// Returns true if the projection was performed.
    ///
    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// Mark the projection as done.
    ///
    /// OCCT: `Done()`.
    pub fn done(&mut self) {
        self.is_done = true;
    }

    /// Returns the type of the projected curve.
    ///
    /// OCCT: `GetType()`.
    pub fn get_type(&self) -> CurveType {
        self.proj_type
    }

    /// Set the curve type.
    ///
    /// OCCT: `SetType(Type)`.
    pub fn set_type(&mut self, typ: CurveType) {
        self.proj_type = typ;
    }

    /// Returns true if the result is periodic.
    ///
    /// OCCT: `IsPeriodic()`.
    pub fn is_periodic(&self) -> bool {
        self.is_periodic
    }

    /// Set periodic flag.
    ///
    /// OCCT: `SetPeriodic()`.
    pub fn set_periodic(&mut self) {
        self.is_periodic = true;
    }

    /// Returns the line result.
    ///
    /// OCCT: `Line()`.
    pub fn line(&self) -> &Line3 {
        &self.lin
    }

    /// Returns the circle result.
    ///
    /// OCCT: `Circle()`.
    pub fn circle(&self) -> &Circle3 {
        &self.circ
    }

    /// Returns the ellipse result.
    ///
    /// OCCT: `Ellipse()`.
    pub fn ellipse(&self) -> &Ellipse3 {
        &self.elips
    }

    /// Returns the hyperbola result.
    ///
    /// OCCT: `Hyperbola()`.
    pub fn hyperbola(&self) -> &Hyperbola3 {
        &self.hypr
    }

    /// Returns the parabola result.
    ///
    /// OCCT: `Parabola()`.
    pub fn parabola(&self) -> &Parabola3 {
        &self.parab
    }

    /// Returns the BSpline/Bezier result.
    ///
    /// OCCT: `BSpline()` / `Bezier()`.
    pub fn bspline(&self) -> Option<&Curve2d> {
        self.bspline.as_ref()
    }

    /// Store a line result.
    pub fn set_line(&mut self, line: Line3) {
        self.proj_type = CurveType::Line;
        self.lin = line;
    }

    /// Store a circle result.
    pub fn set_circle(&mut self, circ: Circle3) {
        self.proj_type = CurveType::Circle;
        self.circ = circ;
    }

    /// Store an ellipse result.
    pub fn set_ellipse(&mut self, elips: Ellipse3) {
        self.proj_type = CurveType::Ellipse;
        self.elips = elips;
    }

    /// Store a BSpline result.
    pub fn set_bspline(&mut self, curve: Curve2d) {
        self.proj_type = CurveType::BSpline;
        self.bspline = Some(curve);
    }

    /// Convert the projection result to a `Curve2d`.
    pub fn to_curve2d(&self) -> Curve2d {
        match self.proj_type {
            CurveType::Line => {
                let l = &self.lin;
                Curve2d::Line(crate::geom::Line2d::new(
                    DVec2::new(l.origin.x, l.origin.y),
                    DVec2::new(l.direction.x, l.direction.y),
                ))
            }
            CurveType::Circle => {
                let c = &self.circ;
                // OCCT ProjLib_Plane/Cylinder::Project(gp_Circ): the projected
                // circle keeps its local frame (gp_Ax22d) — the x/y directions
                // of the 3D circle mapped onto the surface. Dropping the frame
                // (Circle2d::new defaults to the identity axes) parameterizes
                // the pcurve from the surface's u axis instead, so an edge
                // whose arc is not centered on that axis traces the COMPLEMENT
                // half of its circle (bfuse_simple/E1: p1's bottom-cap arc
                // pcurves became the opposite semicircles, misclassifying the
                // EF midpoint and dropping the top-piece common part).
                Curve2d::Circle(crate::geom::Circle2d {
                    center: DVec2::new(c.center.x, c.center.y),
                    x_dir: DVec2::new(c.x_dir.x, c.x_dir.y),
                    y_dir: DVec2::new(c.y_dir.x, c.y_dir.y),
                    radius: c.radius,
                })
            }
            CurveType::Ellipse => {
                let e = &self.elips;
                Curve2d::Ellipse(crate::geom::Ellipse2d {
                    center: DVec2::new(e.center.x, e.center.y),
                    major_dir: DVec2::new(e.major_dir.x, e.major_dir.y),
                    major_radius: e.major_radius,
                    minor_radius: e.minor_radius,
                })
            }
            CurveType::BSpline => {
                self.bspline.clone().unwrap_or(Curve2d::Line(crate::geom::Line2d::new(
                    DVec2::ZERO,
                    DVec2::X,
                )))
            }
            _ => Curve2d::Line(crate::geom::Line2d::new(DVec2::ZERO, DVec2::X)),
        }
    }
}

impl Default for Projector {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// ProjLib_Plane
// ============================================================================

/// Project elementary curves onto a plane.
///
/// OCCT: `ProjLib_Plane`.
pub struct PlaneProjector {
    projector: Projector,
    plane: Plane,
}

impl PlaneProjector {
    /// Create an uninitialized plane projector.
    ///
    /// OCCT: `ProjLib_Plane()`.
    pub fn new() -> Self {
        PlaneProjector {
            projector: Projector::new(),
            plane: Plane::new(DVec3::ZERO, DVec3::Z),
        }
    }

    /// Create a plane projector for the given plane.
    ///
    /// OCCT: `ProjLib_Plane(gp_Pln)`.
    pub fn with_plane(plane: &Plane) -> Self {
        PlaneProjector {
            projector: Projector::new(),
            plane: *plane,
        }
    }

    /// Initialize (or re-initialize) the plane.
    ///
    /// OCCT: `Init(gp_Pln)`.
    pub fn init(&mut self, plane: &Plane) {
        self.plane = *plane;
    }

    /// Project a 3D line onto the plane → 2D line.
    ///
    /// OCCT: `Project(gp_Lin)`.
    pub fn project_line(&mut self, line: &Line3) {
        // Map the line's origin and direction to UV coordinates on the plane
        let o = self.to_uv(line.origin);
        let p1 = self.to_uv(line.origin + line.direction);
        let dir = (p1 - o).normalize_or_zero();
        self.projector.set_line(Line3 {
            origin: DVec3::new(o.x, o.y, 0.0),
            direction: DVec3::new(dir.x, dir.y, 0.0),
        });
        self.projector.done();
    }

    /// Project a 3D circle onto the plane → 2D circle (or ellipse).
    ///
    /// OCCT: `Project(gp_Circ)`.
    pub fn project_circle(&mut self, circle: &Circle3) {
        let center_uv = self.to_uv(circle.center);
        // Project the x_dir and y_dir onto the plane
        let x_uv = self.to_uv(circle.center + circle.x_dir * circle.radius);
        let y_uv = self.to_uv(circle.center + circle.y_dir * circle.radius);

        let x_vec = (x_uv - center_uv) / circle.radius;
        let y_vec = (y_uv - center_uv) / circle.radius;

        // Check if the projection preserves the circle or becomes an ellipse
        let x_len = x_vec.length();
        let y_len = y_vec.length();

        if (x_len - y_len).abs() < 1e-10 && x_vec.dot(y_vec).abs() < 1e-10 {
            // Orthogonal with equal scale → circle
            self.projector.set_circle(Circle3 {
                center: DVec3::new(center_uv.x, center_uv.y, 0.0),
                normal: DVec3::Z,
                x_dir: DVec3::new(x_vec.x / x_len, x_vec.y / x_len, 0.0),
                y_dir: DVec3::new(y_vec.x / y_len, y_vec.y / y_len, 0.0),
                radius: circle.radius * x_len,
            });
        } else {
            // General case → ellipse
            self.projector.set_ellipse(Ellipse3 {
                center: DVec3::new(center_uv.x, center_uv.y, 0.0),
                normal: DVec3::Z,
                major_dir: DVec3::new(x_vec.x, x_vec.y, 0.0).normalize_or_zero(),
                major_radius: circle.radius * x_len,
                minor_radius: circle.radius * y_len,
            });
        }
        self.projector.done();
    }

    /// Project a 3D ellipse onto the plane → 2D ellipse.
    ///
    /// OCCT: `Project(gp_Elips)`.
    pub fn project_ellipse(&mut self, ellipse: &Ellipse3) {
        let center_uv = self.to_uv(ellipse.center);
        let major_uv = self.to_uv(ellipse.center + ellipse.major_dir * ellipse.major_radius);
        let minor_uv = self.to_uv(
            ellipse.center + ellipse.normal.cross(ellipse.major_dir).normalize_or_zero() * ellipse.minor_radius,
        );

        let major_vec = (major_uv - center_uv) / ellipse.major_radius;
        let minor_vec = (minor_uv - center_uv) / ellipse.minor_radius;

        self.projector.set_ellipse(Ellipse3 {
            center: DVec3::new(center_uv.x, center_uv.y, 0.0),
            normal: DVec3::Z,
            major_dir: DVec3::new(major_vec.x, major_vec.y, 0.0).normalize_or_zero(),
            major_radius: ellipse.major_radius * major_vec.length(),
            minor_radius: ellipse.minor_radius * minor_vec.length(),
        });
        self.projector.done();
    }

    /// Project a 3D parabola onto the plane → 2D parabola.
    ///
    /// OCCT: `Project(gp_Parab)`.
    pub fn project_parabola(&mut self, parabola: &Parabola3) {
        let vertex_uv = self.to_uv(parabola.vertex);
        let axis_uv = self.to_uv(parabola.vertex + parabola.axis_dir);

        let axis_dir_2d = (axis_uv - vertex_uv).normalize_or_zero();

        self.projector.set_type(CurveType::Parabola);
        self.projector.parab = Parabola3 {
            vertex: DVec3::new(vertex_uv.x, vertex_uv.y, 0.0),
            normal: DVec3::Z,
            axis_dir: DVec3::new(axis_dir_2d.x, axis_dir_2d.y, 0.0),
            focal_param: parabola.focal_param * axis_dir_2d.length(),
        };
        self.projector.done();
    }

    /// Project a 3D hyperbola onto the plane → 2D hyperbola.
    ///
    /// OCCT: `Project(gp_Hypr)`.
    pub fn project_hyperbola(&mut self, hyperbola: &Hyperbola3) {
        let center_uv = self.to_uv(hyperbola.center);
        let major_uv = self.to_uv(hyperbola.center + hyperbola.major_dir);

        let major_dir_2d = (major_uv - center_uv).normalize_or_zero();
        let scale = major_dir_2d.length();

        self.projector.set_type(CurveType::Hyperbola);
        self.projector.hypr = Hyperbola3 {
            center: DVec3::new(center_uv.x, center_uv.y, 0.0),
            normal: DVec3::Z,
            major_dir: DVec3::new(major_dir_2d.x, major_dir_2d.y, 0.0),
            semi_major: hyperbola.semi_major * scale,
            semi_minor: hyperbola.semi_minor * scale,
        };
        self.projector.done();
    }

    /// Returns the projector (result).
    pub fn projector(&self) -> &Projector {
        &self.projector
    }

    fn to_uv(&self, p: DVec3) -> DVec2 {
        let d = p - self.plane.origin;
        DVec2::new(d.dot(self.plane.u_dir), d.dot(self.plane.v_dir))
    }
}

impl Default for PlaneProjector {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// ProjLib_Cylinder
// ============================================================================

/// Project elementary curves onto a cylinder (develop into UV rectangle).
///
/// OCCT: `ProjLib_Cylinder`.
pub struct CylinderProjector {
    projector: Projector,
    cylinder: CylindricalSurface,
}

impl CylinderProjector {
    pub fn new() -> Self {
        CylinderProjector {
            projector: Projector::new(),
            cylinder: CylindricalSurface {
                origin: DVec3::ZERO,
                axis: DVec3::Z,
                radius: 1.0,
                ref_dir: DVec3::X,
            },
        }
    }

    pub fn with_cylinder(cylinder: &CylindricalSurface) -> Self {
        CylinderProjector {
            projector: Projector::new(),
            cylinder: *cylinder,
        }
    }

    pub fn init(&mut self, cylinder: &CylindricalSurface) {
        self.cylinder = *cylinder;
    }

    /// Project a 3D line onto the cylinder surface parameter space.
    ///
    /// OCCT: `Project(gp_Lin)`.
    pub fn project_line(&mut self, line: &Line3) {
        let o = self.to_uv(line.origin);
        let p1 = self.to_uv(line.origin + line.direction);
        let dir = (p1 - o).normalize_or_zero();
        self.projector.set_line(Line3 {
            origin: DVec3::new(o.x, o.y, 0.0),
            direction: DVec3::new(dir.x, dir.y, 0.0),
        });
        self.projector.done();
        if dir.x.abs() < 1e-12 {
            self.projector.set_periodic();
        }
    }

    /// Project a 3D circle onto the cylinder parameter space.
    ///
    /// OCCT: `Project(gp_Circ)`.
    pub fn project_circle(&mut self, circle: &Circle3) {
        let center_uv = self.to_uv(circle.center);
        let r = circle.radius / self.cylinder.radius;

        // Check if the circle lies in a plane parallel to the cylinder axis
        let cos_angle = circle.normal.dot(self.cylinder.axis).abs();
        if cos_angle < 1e-12 {
            // Circle in plane parallel to axis → ellipse in UV
            self.projector.set_ellipse(Ellipse3 {
                center: DVec3::new(center_uv.x, center_uv.y, 0.0),
                normal: DVec3::Z,
                major_dir: DVec3::X,
                major_radius: r,
                minor_radius: r,
            });
        } else {
            // General case: approximate with BSpline
            let n = 16;
            let mut pts = Vec::with_capacity(n);
            for i in 0..n {
                let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
                let p = circle.center
                    + circle.x_dir * circle.radius * theta.cos()
                    + circle.y_dir * circle.radius * theta.sin();
                pts.push(self.to_uv(p));
            }
            self.projector.set_bspline(Curve2d::BSpline(
                crate::geom::BSplineCurve2::approximate_closed(&pts),
            ));
        }
        self.projector.done();
    }

    /// Project a 3D ellipse onto the cylinder parameter space.
    ///
    /// OCCT: `Project(gp_Elips)`.
    pub fn project_ellipse(&mut self, ellipse: &Ellipse3) {
        let n = 32;
        let mut pts = Vec::with_capacity(n);
        for i in 0..n {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            let p = ellipse.center
                + ellipse.major_dir * ellipse.major_radius * theta.cos()
                + ellipse.normal.cross(ellipse.major_dir).normalize_or_zero() * ellipse.minor_radius * theta.sin();
            pts.push(self.to_uv(p));
        }
        self.projector.set_bspline(Curve2d::BSpline(
            crate::geom::BSplineCurve2::approximate_closed(&pts),
        ));
        self.projector.done();
    }

    /// Returns the projector (result).
    pub fn projector(&self) -> &Projector {
        &self.projector
    }

    fn to_uv(&self, p: DVec3) -> DVec2 {
        let d = p - self.cylinder.origin;
        let along = d.dot(self.cylinder.axis);
        let radial = d - self.cylinder.axis * along;
        let u = radial.dot(self.cylinder.ref_dir).atan2(
            self.cylinder.axis.cross(self.cylinder.ref_dir).normalize_or_zero().dot(radial),
        );
        DVec2::new(u, along)
    }
}

impl Default for CylinderProjector {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// ProjLib_Sphere
// ============================================================================

/// Project elementary curves onto a sphere (into UV parameter space).
///
/// OCCT: `ProjLib_Sphere`.
pub struct SphereProjector {
    projector: Projector,
    sphere: SphericalSurface,
}

impl SphereProjector {
    pub fn new() -> Self {
        SphereProjector {
            projector: Projector::new(),
            sphere: SphericalSurface {
                center: DVec3::ZERO,
                axis: DVec3::Z,
                radius: 1.0,
                ref_dir: DVec3::X,
            },
        }
    }

    pub fn with_sphere(sphere: &SphericalSurface) -> Self {
        SphereProjector {
            projector: Projector::new(),
            sphere: *sphere,
        }
    }

    fn to_uv(&self, p: DVec3) -> DVec2 {
        let d = p - self.sphere.center;
        let r = d.length();
        if r < 1e-15 {
            return DVec2::ZERO;
        }
        let d_n = d / r;
        let u = d_n.dot(self.sphere.ref_dir_perp()).atan2(d_n.dot(self.sphere.ref_dir));
        let v = d_n.dot(self.sphere.axis).acos();
        DVec2::new(u, v)
    }

    fn project_generic(&mut self, n: usize, points: &[DVec3]) {
        let mut pts: Vec<DVec2> = points.iter().map(|&p| self.to_uv(p)).collect();
        if pts.is_empty() {
            return;
        }
        // Detect periodicity in U
        let mut u_min = pts[0].x;
        let mut u_max = pts[0].x;
        for pt in &pts {
            u_min = u_min.min(pt.x);
            u_max = u_max.max(pt.x);
        }
        if u_max - u_min > std::f64::consts::PI {
            // Unwrap U coordinates
            for pt in &mut pts {
                if pt.x < 0.0 {
                    pt.x += 2.0 * std::f64::consts::PI;
                }
            }
        }
        self.projector.set_bspline(Curve2d::BSpline(
            crate::geom::BSplineCurve2::approximate(&pts),
        ));
        self.projector.done();
    }

    pub fn project_line(&mut self, line: &Line3) {
        let n = 32;
        let pts: Vec<DVec3> = (0..n).map(|i| {
            let t = (i as f64) / (n as f64) * 4.0;
            line.origin + line.direction * t
        }).collect();
        self.project_generic(n, &pts);
    }

    pub fn project_circle(&mut self, circle: &Circle3) {
        let n = 32;
        let pts: Vec<DVec3> = (0..n).map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            circle.center + circle.x_dir * circle.radius * theta.cos()
                + circle.y_dir * circle.radius * theta.sin()
        }).collect();
        self.project_generic(n, &pts);
    }

    pub fn projector(&self) -> &Projector {
        &self.projector
    }
}

// ============================================================================
// ProjLib_Cone
// ============================================================================

/// Project elementary curves onto a cone (develop into UV space).
///
/// OCCT: `ProjLib_Cone`.
pub struct ConeProjector {
    projector: Projector,
    cone: ConicalSurface,
}

impl ConeProjector {
    pub fn new() -> Self {
        ConeProjector {
            projector: Projector::new(),
            cone: ConicalSurface::new(DVec3::ZERO, DVec3::Z, 1.0, 0.25),
        }
    }

    pub fn with_cone(cone: &ConicalSurface) -> Self {
        ConeProjector {
            projector: Projector::new(),
            cone: *cone,
        }
    }

    fn to_uv(&self, p: DVec3) -> DVec2 {
        self.cone.world_to_uv(p)
    }

    fn project_generic(&mut self, n: usize, points: &[DVec3]) {
        let pts: Vec<DVec2> = points.iter().map(|&p| self.to_uv(p)).collect();
        self.projector.set_bspline(Curve2d::BSpline(
            crate::geom::BSplineCurve2::approximate(&pts),
        ));
        self.projector.done();
    }

    pub fn project_line(&mut self, line: &Line3) {
        let n = 32;
        let pts: Vec<DVec3> = (0..n).map(|i| {
            let t = (i as f64) / (n as f64) * 4.0;
            line.origin + line.direction * t
        }).collect();
        self.project_generic(n, &pts);
    }

    pub fn project_circle(&mut self, circle: &Circle3) {
        let n = 32;
        let pts: Vec<DVec3> = (0..n).map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            circle.center + circle.x_dir * circle.radius * theta.cos()
                + circle.y_dir * circle.radius * theta.sin()
        }).collect();
        self.project_generic(n, &pts);
    }

    pub fn projector(&self) -> &Projector {
        &self.projector
    }
}

// ============================================================================
// ProjLib_Torus
// ============================================================================

/// Project curves onto a torus (into UV parameter space).
///
/// OCCT: `ProjLib_Torus`.
pub struct TorusProjector {
    projector: Projector,
    torus: ToroidalSurface,
}

impl TorusProjector {
    pub fn new() -> Self {
        TorusProjector {
            projector: Projector::new(),
            torus: ToroidalSurface {
                center: DVec3::ZERO,
                axis: DVec3::Z,
                major_radius: 3.0,
                minor_radius: 1.0,
            },
        }
    }

    pub fn with_torus(torus: &ToroidalSurface) -> Self {
        TorusProjector {
            projector: Projector::new(),
            torus: *torus,
        }
    }

    fn to_uv(&self, p: DVec3) -> DVec2 {
        let d = p - self.torus.center;
        let along = d.dot(self.torus.axis);
        let radial = d - self.torus.axis * along;
        let r = radial.length();
        let u = radial.y.atan2(radial.x);
        let v = if r < 1e-15 {
            0.0
        } else {
            let local_r = r - self.torus.major_radius;
            let v_angle = local_r.atan2(along);
            v_angle
        };
        DVec2::new(u, v)
    }

    fn project_generic(&mut self, n: usize, points: &[DVec3]) {
        let pts: Vec<DVec2> = points.iter().map(|&p| self.to_uv(p)).collect();
        self.projector.set_bspline(Curve2d::BSpline(
            crate::geom::BSplineCurve2::approximate(&pts),
        ));
        self.projector.done();
    }

    pub fn project_line(&mut self, line: &Line3) {
        let n = 32;
        let pts: Vec<DVec3> = (0..n).map(|i| {
            let t = (i as f64) / (n as f64) * 4.0;
            line.origin + line.direction * t
        }).collect();
        self.project_generic(n, &pts);
    }

    pub fn project_circle(&mut self, circle: &Circle3) {
        let n = 32;
        let pts: Vec<DVec3> = (0..n).map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            circle.center + circle.x_dir * circle.radius * theta.cos()
                + circle.y_dir * circle.radius * theta.sin()
        }).collect();
        self.project_generic(n, &pts);
    }

    pub fn projector(&self) -> &Projector {
        &self.projector
    }
}

// ============================================================================
// Free-function convenience API (Rust-style, complements OCCT-aligned classes)
// ============================================================================

/// Project a 3D curve onto a surface, returning the 2D pcurve.
pub fn project_on_surface(curve: &Curve3, surface: &Surface3) -> Option<Curve2d> {
    match (curve, surface) {
        (Curve3::Line(l), Surface3::Plane(p)) => {
            let mut proj = PlaneProjector::with_plane(p);
            proj.project_line(l);
            Some(proj.projector().to_curve2d())
        }
        (Curve3::Line(l), Surface3::Cylinder(c)) => {
            let mut proj = CylinderProjector::with_cylinder(c);
            proj.project_line(l);
            Some(proj.projector().to_curve2d())
        }
        (Curve3::Circle(c), Surface3::Plane(p)) => {
            let mut proj = PlaneProjector::with_plane(p);
            proj.project_circle(c);
            Some(proj.projector().to_curve2d())
        }
        (Curve3::Circle(c), Surface3::Cylinder(cyl)) => {
            let mut proj = CylinderProjector::with_cylinder(cyl);
            proj.project_circle(c);
            Some(proj.projector().to_curve2d())
        }
        _ => None, // Complex cases require numerical approximation
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Curve2dEval;

    #[test]
    fn test_plane_projector_line() {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let line = Line3::new(DVec3::new(1.0, 2.0, 3.0), DVec3::X);
        let mut proj = PlaneProjector::with_plane(&plane);
        proj.project_line(&line);
        assert!(proj.projector().is_done());
        assert_eq!(proj.projector().get_type(), CurveType::Line);
    }

    #[test]
    fn test_plane_projector_circle() {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let circle = Circle3::new(DVec3::new(0.0, 0.0, 5.0), DVec3::Z, 3.0);
        let mut proj = PlaneProjector::with_plane(&plane);
        proj.project_circle(&circle);
        assert!(proj.projector().is_done());
        assert_eq!(proj.projector().get_type(), CurveType::Circle);
    }

    #[test]
    fn test_cylinder_projector_line() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        let line = Line3::new(DVec3::new(5.0, 0.0, 0.0), DVec3::Z);
        let mut proj = CylinderProjector::with_cylinder(&cyl);
        proj.project_line(&line);
        assert!(proj.projector().is_done());
    }

    #[test]
    fn test_project_on_surface_line_plane() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let result = project_on_surface(&line, &plane);
        assert!(result.is_some());
        let pcurve = result.unwrap();
        let p0 = pcurve.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-10);
    }
}
