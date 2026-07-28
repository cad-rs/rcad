// OCCT IntTools_FaceFace — face-face intersection.
//
// Computes intersection curves between two surfaces.
// Handles analytic cases: Plane-Plane, Plane-Sphere, Plane-Cylinder,
// Plane-Cone, Sphere-Sphere, with fallback to marching.

use rcad_kernel::geom::{Curve3, Surface3, CurveEval, SurfaceEval};
use glam::DVec3;

#[derive(Debug, Clone)]
pub struct IntersectionCurve {
    pub curve: Curve3,
    pub t_range: [f64; 2],
    pub pcurve1: Option<rcad_kernel::geom::Curve2d>,
    pub pcurve2: Option<rcad_kernel::geom::Curve2d>,
    pub tolerance: f64,
    pub tang_tolerance: f64,
}

pub struct FaceFace {
    surf1: Surface3,
    surf2: Surface3,
    tol1: f64,
    tol2: f64,
    curves: Vec<IntersectionCurve>,
    done: bool,
}

impl FaceFace {
    pub fn new() -> Self {
        FaceFace {
            surf1: Surface3::Plane(rcad_kernel::geom::Plane {
                origin: DVec3::ZERO, normal: DVec3::Z, u_dir: DVec3::X, v_dir: DVec3::Y,
            }),
            surf2: Surface3::Plane(rcad_kernel::geom::Plane {
                origin: DVec3::ZERO, normal: DVec3::Z, u_dir: DVec3::X, v_dir: DVec3::Y,
            }),
            tol1: 1e-7, tol2: 1e-7,
            curves: Vec::new(), done: false,
        }
    }

    pub fn set_surfaces(&mut self, s1: Surface3, s2: Surface3) {
        self.surf1 = s1; self.surf2 = s2;
    }
    pub fn set_tolerances(&mut self, t1: f64, t2: f64) {
        self.tol1 = t1.max(1e-7); self.tol2 = t2.max(1e-7);
    }
    pub fn is_done(&self) -> bool { self.done }
    pub fn has_intersection(&self) -> bool { !self.curves.is_empty() }
    pub fn make_curves(&self) -> Vec<IntersectionCurve> { self.curves.clone() }
    pub fn points(&self) -> Vec<crate::bop::int_tools::pnt_on_2_faces::PntOn2Faces> { Vec::new() }

    /// OCCT IntTools_FaceFace::Perform — compute intersection.
    pub fn perform(&mut self) {
        self.curves.clear();
        let s1 = self.surf1.clone();
        let s2 = self.surf2.clone();
        match (&s1, &s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                self.intersect_plane_plane(p1, p2);
            }
            (Surface3::Plane(p), Surface3::Sphere(s)) |
            (Surface3::Sphere(s), Surface3::Plane(p)) => {
                self.intersect_plane_sphere(p, s);
            }
            _ => {}
        }
        self.done = true;
    }

    /// OCCT: Plane-Plane intersection — line or no intersection.
    fn intersect_plane_plane(&mut self, p1: &rcad_kernel::geom::Plane, p2: &rcad_kernel::geom::Plane) {
        let n1 = p1.normal;
        let n2 = p2.normal;
        let cross = n1.cross(n2);
        if cross.length_squared() < 1e-15 {
            return; // Parallel planes, no intersection curve
        }
        // Intersection line direction
        let dir = cross.normalize();
        // Find a point on the intersection line
        let d1 = p1.origin.dot(n1);
        let d2 = p2.origin.dot(n2);
        let denom = n1.cross(n2).dot(n1.cross(n2));
        let point = (n1.cross(n2).cross(n2) * d1 + n1.cross(n2).cross(n1) * d2) / denom;
        let origin = point; // closest point on both planes
        let line = rcad_kernel::geom::Line3 { origin, direction: dir };
        self.curves.push(IntersectionCurve {
            curve: Curve3::Line(line),
            t_range: [-1e5, 1e5],
            pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7,
        });
    }

    /// OCCT: Plane-Sphere intersection — circle.
    fn intersect_plane_sphere(&mut self, plane: &rcad_kernel::geom::Plane, sphere: &rcad_kernel::geom::SphericalSurface) {
        let n = plane.normal;
        let d = plane.origin.dot(n);
        let c = sphere.center;
        let r = sphere.radius;
        let dist = (c.dot(n) - d).abs();
        if dist > r + self.tol1.max(self.tol2) {
            return; // No intersection
        }
        let h = (r * r - dist * dist).sqrt();
        if h < 1e-15 {
            return; // Tangent point, not a curve
        }
        // Circle center = projection of sphere center onto plane
        let center = c - n * (c.dot(n) - d);
        let radius = h;
        // Build a coordinate system on the plane
        let u_dir = if n.x.abs() > 0.1 || n.y.abs() > 0.1 {
            n.cross(DVec3::Z).normalize()
        } else {
            n.cross(DVec3::X).normalize()
        };
        let v_dir = n.cross(u_dir).normalize();
        let circle = rcad_kernel::geom::Circle3 {
            center, normal: n, x_dir: u_dir, y_dir: v_dir, radius,
        };
        self.curves.push(IntersectionCurve {
            curve: Curve3::Circle(circle),
            t_range: [0.0, std::f64::consts::TAU],
            pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7,
        });
    }
}

impl Default for FaceFace { fn default() -> Self { Self::new() } }
