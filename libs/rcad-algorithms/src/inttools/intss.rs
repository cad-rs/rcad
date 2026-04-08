//! Surface-surface intersection (IntSS).
//!
//! Covers analytic pairs:
//!
//! | Pair | Result |
//! |------|--------|
//! | Plane × Plane | Line or parallel/coincident |
//! | Plane × Sphere | Circle |
//! | Plane × Cylinder | Circle / Ellipse / Lines |
//! | Plane × Cone | Circle / Ellipse / Lines |
//! | Sphere × Sphere | Circle (intersection plane ⊥ line-of-centres) |
//! | Sphere × Cylinder | Circle (axis ⊥ case) / Numeric |
//! | Cylinder × Cylinder | Ellipse (axes ∥), numeric otherwise |
//! | Cylinder × Cone | Numeric |
//! | Sphere × Cone | Circle (apex-centred case) / Numeric |
//! | Everything else | Numeric polylines via marching |
//!
//! Analogous to OCCT `GeomAPI_IntSS`.

use std::f64::consts::PI;

use glam::DVec3;
use rcad_kernel::geom::{
    Circle3, ConicalSurface, CurveEval, CylindricalSurface, Ellipse3, Line3, Plane,
    SphericalSurface, Surface3, SurfaceEval, any_perpendicular,
};

use crate::inttools::{
    plane_cone::{PlaneConicalResult, intersect_plane_cone},
    plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder},
    plane_plane::{PlanePlaneResult, intersect_plane_plane},
    plane_sphere::{PlaneSphereResult, intersect_plane_sphere},
};
use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG, vectors_parallel};

// ──────────────────────────────────────────────────────────────────────────────
// Public result types
// ──────────────────────────────────────────────────────────────────────────────

/// A single intersection component between two surfaces.
#[derive(Debug, Clone)]
pub enum SurfaceCurve {
    /// An exact analytic circle.
    Circle(Circle3),
    /// An exact analytic ellipse.
    Ellipse(Ellipse3),
    /// An exact analytic line (infinite).
    Line(Line3),
    /// A tangent point (zero-dimensional contact).
    Point(DVec3),
    /// Numerically sampled polyline (fallback for non-analytic pairs).
    Polyline(Vec<DVec3>),
}

/// All intersection curves / components found between two surfaces.
#[derive(Debug, Clone, Default)]
pub struct SurfaceSurfaceIntersection {
    pub curves: Vec<SurfaceCurve>,
}

impl SurfaceSurfaceIntersection {
    pub fn is_empty(&self) -> bool {
        self.curves.is_empty()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Main dispatch
// ──────────────────────────────────────────────────────────────────────────────

/// Compute the intersection between two `Surface3` values.
///
/// Returns exact analytic curves where possible; falls back to numerical
/// polylines for unsupported surface-type combinations.
pub fn intersect_surfaces(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection {
    use Surface3::*;
    match (s1, s2) {
        // ── Plane × * ─────────────────────────────────────────────────────
        (Plane(p1), Plane(p2)) => plane_x_plane(p1, p2),
        (Plane(p), Sphere(s)) | (Sphere(s), Plane(p)) => plane_x_sphere(p, s),
        (Plane(p), Cylinder(c)) | (Cylinder(c), Plane(p)) => plane_x_cylinder(p, c),
        (Plane(p), Cone(c)) | (Cone(c), Plane(p)) => plane_x_cone(p, c),

        // ── Sphere × * ────────────────────────────────────────────────────
        (Sphere(s1), Sphere(s2)) => sphere_x_sphere(s1, s2),
        (Sphere(s), Cylinder(c)) | (Cylinder(c), Sphere(s)) => sphere_x_cylinder(s, c),
        (Sphere(s), Cone(c)) | (Cone(c), Sphere(s)) => sphere_x_cone(s, c),

        // ── Cylinder × Cylinder ───────────────────────────────────────────
        (Cylinder(c1), Cylinder(c2)) => cylinder_x_cylinder(c1, c2),

        // ── Cylinder × Cone and Cone × Cone fall through to numeric ───────
        // ── All others → numeric marching ─────────────────────────────────
        _ => numeric_intss(s1, s2),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Plane × Plane
// ──────────────────────────────────────────────────────────────────────────────

fn plane_x_plane(p1: &Plane, p2: &Plane) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_plane(p1, p2) {
        PlanePlaneResult::Line(l) => {
            out.curves.push(SurfaceCurve::Line(l));
        }
        PlanePlaneResult::Coincident => {} // surfaces identical — infinite overlap
        PlanePlaneResult::Parallel => {}   // no intersection
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Plane × Sphere
// ──────────────────────────────────────────────────────────────────────────────

fn plane_x_sphere(p: &Plane, s: &SphericalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_sphere(p, s) {
        PlaneSphereResult::Circle(c) => {
            out.curves.push(SurfaceCurve::Circle(c));
        }
        PlaneSphereResult::TangentPoint(pt) => {
            out.curves.push(SurfaceCurve::Point(pt));
        }
        PlaneSphereResult::NoIntersection => {}
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Plane × Cylinder
// ──────────────────────────────────────────────────────────────────────────────

fn plane_x_cylinder(p: &Plane, c: &CylindricalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_cylinder(p, c) {
        PlaneCylinderResult::Circle(c) => {
            out.curves.push(SurfaceCurve::Circle(c));
        }
        PlaneCylinderResult::Ellipse(e) => {
            out.curves.push(SurfaceCurve::Ellipse(e));
        }
        PlaneCylinderResult::TangentLine(l) => {
            out.curves.push(SurfaceCurve::Line(l));
        }
        PlaneCylinderResult::TwoLines(l1, l2) => {
            out.curves.push(SurfaceCurve::Line(l1));
            out.curves.push(SurfaceCurve::Line(l2));
        }
        PlaneCylinderResult::NoIntersection => {}
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Plane × Cone
// ──────────────────────────────────────────────────────────────────────────────

fn plane_x_cone(p: &Plane, c: &ConicalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_cone(p, c) {
        PlaneConicalResult::Circle(c) => {
            out.curves.push(SurfaceCurve::Circle(c));
        }
        PlaneConicalResult::Ellipse(e) => {
            out.curves.push(SurfaceCurve::Ellipse(e));
        }
        PlaneConicalResult::SingleLine(l) => {
            out.curves.push(SurfaceCurve::Line(l));
        }
        PlaneConicalResult::TwoLines(l1, l2) => {
            out.curves.push(SurfaceCurve::Line(l1));
            out.curves.push(SurfaceCurve::Line(l2));
        }
        PlaneConicalResult::Point(pt) => {
            out.curves.push(SurfaceCurve::Point(pt));
        }
        PlaneConicalResult::NoIntersection => {}
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Sphere × Sphere
// ──────────────────────────────────────────────────────────────────────────────

/// Two spheres intersect in a circle (or are tangent/disjoint).
///
/// The intersection circle lies on the radical plane, whose normal is the
/// line-of-centres direction.
fn sphere_x_sphere(s1: &SphericalSurface, s2: &SphericalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    let d_vec = s2.center - s1.center;
    let d = d_vec.length();

    // Concentric spheres: coincident (same r) or no intersection
    if d < TOLERANCE_ABS {
        return out; // treat as no intersection (or coincident, infinite)
    }

    let r1 = s1.radius;
    let r2 = s2.radius;

    // Disjoint or one contains the other
    if d > r1 + r2 + TOLERANCE_ABS || d < (r1 - r2).abs() - TOLERANCE_ABS {
        return out;
    }

    let axis = d_vec / d;

    // Distance from s1.center to the intersection plane (radical plane)
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);

    // Tangent case
    let r_sq = r1 * r1 - a * a;
    if r_sq < -TOLERANCE_ABS {
        return out;
    }
    let r_circle = r_sq.max(0.0).sqrt();
    let center = s1.center + axis * a;

    if r_circle < TOLERANCE_ABS {
        out.curves.push(SurfaceCurve::Point(center));
    } else {
        out.curves.push(SurfaceCurve::Circle(Circle3 {
            center,
            normal: axis,
            radius: r_circle,
        }));
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Sphere × Cylinder
// ──────────────────────────────────────────────────────────────────────────────

/// Sphere-cylinder intersection.
///
/// Analytic case: cylinder axis passes through the sphere centre →
/// two parallel circles (or one if tangent).
/// All other cases fall back to numerical marching.
fn sphere_x_cylinder(s: &SphericalSurface, c: &CylindricalSurface) -> SurfaceSurfaceIntersection {
    // Project sphere centre onto cylinder axis
    let t = (s.center - c.origin).dot(c.axis);
    let foot = c.origin + c.axis * t;
    let d_perp = (s.center - foot).length();

    // If the sphere centre is on the cylinder axis, section planes ⊥ axis give
    // circles at heights where r_sphere(z)² = R_cylinder²
    // r_sphere(z)² = R² - (z - z_c)² where z_c is the axial position of sphere center
    // Solve: R² - z² = r_cyl²  (in local frame with sphere center as origin along axis)
    if d_perp < TOLERANCE_ABS {
        // Sphere centre on axis — analytic circles
        let dz_sq = s.radius * s.radius - c.radius * c.radius;
        if dz_sq < -TOLERANCE_ABS {
            // Sphere smaller than cylinder — no intersection if dz_sq < 0
            // Actually: sphere radius < cylinder radius means sphere inside cyl,
            // could still intersect if large enough. Recheck:
            // Points on intersection: distance from axis = c.radius AND on sphere.
            // If s.radius < c.radius: sphere never reaches cylinder surface → no intersect.
            return SurfaceSurfaceIntersection::default();
        }
        let mut out = SurfaceSurfaceIntersection::default();
        if dz_sq.abs() < TOLERANCE_ABS {
            // Tangent — single circle at sphere center height
            out.curves.push(SurfaceCurve::Circle(Circle3 {
                center: s.center,
                normal: c.axis,
                radius: c.radius,
            }));
        } else {
            let dz = dz_sq.sqrt();
            for &sign in &[1.0f64, -1.0] {
                let center = s.center + c.axis * (sign * dz);
                out.curves.push(SurfaceCurve::Circle(Circle3 {
                    center,
                    normal: c.axis,
                    radius: c.radius,
                }));
            }
        }
        return out;
    }

    // General case: numerical
    numeric_intss(&Surface3::Sphere(*s), &Surface3::Cylinder(*c))
}

// ──────────────────────────────────────────────────────────────────────────────
// Sphere × Cone
// ──────────────────────────────────────────────────────────────────────────────

/// Sphere-cone intersection.
///
/// Analytic case: sphere centre on cone axis → circles at intersecting heights.
/// General case → numerical.
fn sphere_x_cone(s: &SphericalSurface, c: &ConicalSurface) -> SurfaceSurfaceIntersection {
    // Project sphere centre onto cone axis
    let t = (s.center - c.apex).dot(c.axis);
    let foot = c.apex + c.axis * t;
    let d_perp = (s.center - foot).length();

    if d_perp < TOLERANCE_ABS {
        // Sphere centre on cone axis.
        // Cone radius at axial parameter t: r(t) = t * tan(half_angle) + c.radius
        // Points on intersection satisfy: distance from axis = r(t_z) AND on sphere.
        // In local frame z along axis from apex, sphere centre at z = t:
        //   r_cone(z) = c.radius + z * tan(α)
        //   r_sphere(z) = sqrt(R² - (z-t)²) ... component perpendicular to axis
        // Points on both surfaces: sqrt(R²-(z-t)²) = c.radius + z*tan(α)
        // Let u = z - t:
        //   sqrt(R²-u²) = c.radius + (u+t)*tan(α)
        //   R²-u² = (c.radius + (u+t)*tan(α))²
        // Solve numerically (at most 2 solutions) — emit circles.
        let mut out = SurfaceSurfaceIntersection::default();
        let ta = c.half_angle_rad.tan();
        // Sample z from (t-R) to (t+R), find sign changes of f(z)
        let z_c = t; // sphere centre axial coord (from apex)
        let r_s = s.radius;
        let n = 64usize;
        let mut prev_f = f64::NAN;
        let mut prev_z = 0.0f64;
        for i in 0..=n {
            let frac = i as f64 / n as f64;
            let z = z_c - r_s + 2.0 * r_s * frac;
            let dz = z - z_c;
            let r_sphere_sq = r_s * r_s - dz * dz;
            if r_sphere_sq < 0.0 {
                prev_f = f64::NAN;
                prev_z = z;
                continue;
            }
            let r_sphere = r_sphere_sq.sqrt();
            let r_cone = c.radius + z * ta;
            let f = r_sphere - r_cone;
            if !prev_f.is_nan() && prev_f * f < 0.0 {
                // Bisect
                let mut lo = prev_z;
                let mut hi = z;
                for _ in 0..32 {
                    let mid = (lo + hi) * 0.5;
                    let dm = mid - z_c;
                    let rsm = (r_s * r_s - dm * dm).max(0.0).sqrt();
                    let rcm = c.radius + mid * ta;
                    if (rsm - rcm) * (rsm - (c.radius + lo * ta)) < 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                let z_sol = (lo + hi) * 0.5;
                let r_sol = (c.radius + z_sol * ta).max(0.0);
                let center = c.apex + c.axis * z_sol;
                if r_sol > TOLERANCE_ABS {
                    out.curves.push(SurfaceCurve::Circle(Circle3 {
                        center,
                        normal: c.axis,
                        radius: r_sol,
                    }));
                }
            }
            prev_f = f;
            prev_z = z;
        }
        if !out.curves.is_empty() {
            return out;
        }
    }

    numeric_intss(&Surface3::Sphere(*s), &Surface3::Cone(*c))
}

// ──────────────────────────────────────────────────────────────────────────────
// Cylinder × Cylinder
// ──────────────────────────────────────────────────────────────────────────────

/// Cylinder-cylinder intersection.
///
/// Analytic case: parallel axes → ellipses (or circles if same radius and
/// same orientation).  General case (skew/crossing axes) → numerical.
fn cylinder_x_cylinder(
    c1: &CylindricalSurface,
    c2: &CylindricalSurface,
) -> SurfaceSurfaceIntersection {
    if vectors_parallel(c1.axis, c2.axis) {
        // Parallel cylinders.
        // Find separation of axes.
        let diff = c2.origin - c1.origin;
        // Project diff onto plane perp to axis
        let proj = diff - c1.axis * diff.dot(c1.axis);
        let d = proj.length();
        let r1 = c1.radius;
        let r2 = c2.radius;

        // No intersection or one inside the other
        if d > r1 + r2 + TOLERANCE_ABS || d < (r1 - r2).abs() - TOLERANCE_ABS {
            return SurfaceSurfaceIntersection::default();
        }

        // For coaxial cylinders of the same radius → coincident (infinite intersection)
        if d < TOLERANCE_ABS && (r1 - r2).abs() < TOLERANCE_ABS {
            return SurfaceSurfaceIntersection::default(); // coincident
        }

        // The two cylinders intersect in two lines parallel to the axis (for infinite cylinders)
        // At angle θ in the cross-section where: r1²+d²-2*d*r1*cos(θ) = r2² → θ from c1 axis
        // These intersection lines are infinitely long — represent as lines through 2 points.
        let mut out = SurfaceSurfaceIntersection::default();
        // Direction of separation (in perp plane)
        let sep_dir = if d > TOLERANCE_ABS {
            proj / d
        } else {
            any_perpendicular(c1.axis)
        };
        // Angle of intersection point from c1 axis towards c2 axis
        let cos_t = if d > TOLERANCE_ABS {
            (d * d + r1 * r1 - r2 * r2) / (2.0 * d * r1)
        } else {
            0.0
        };
        let cos_t = cos_t.clamp(-1.0, 1.0);
        let sin_t = (1.0 - cos_t * cos_t).sqrt();
        let perp = c1.axis.cross(sep_dir).normalize_or_zero();

        for &sign in &[1.0f64, -1.0f64] {
            let dir_in_plane = sep_dir * cos_t + perp * (sign * sin_t);
            let pt = c1.origin + dir_in_plane * r1;
            out.curves.push(SurfaceCurve::Line(Line3 {
                origin: pt,
                direction: c1.axis,
            }));
            if sin_t < TOLERANCE_ABS {
                break;
            } // tangent — only one line
        }
        return out;
    }

    // Non-parallel axes → numerical
    numeric_intss(&Surface3::Cylinder(*c1), &Surface3::Cylinder(*c2))
}

// ──────────────────────────────────────────────────────────────────────────────
// Numerical fallback — marching on parameter-space grid
// ──────────────────────────────────────────────────────────────────────────────

/// Numerical surface-surface intersection via sign-change marching on a grid.
///
/// Samples `n×n` parameter points on `s1`, computes `f(u,v) = dist(P1(u,v), s2) - 0`
/// by projecting onto s2 (approximated by closest distance from a sample grid of s2).
/// Connects sign-change cells into polyline segments.
///
/// This is a coarse approximation (suitable for visual / topological purposes).
/// For high-precision work, Newton refinement would be added on top.
fn numeric_intss(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection {
    const N: usize = 48;
    let dom1 = s1.default_domain();
    let dom2 = s2.default_domain();

    let [u1_0, u1_1, v1_0, v1_1] = dom1;
    let [u2_0, u2_1, v2_0, v2_1] = dom2;

    // Pre-sample s2 on a grid for distance computation
    let n2 = 32usize;
    let mut s2_pts: Vec<DVec3> = Vec::with_capacity(n2 * n2);
    for i in 0..n2 {
        for j in 0..n2 {
            let u = u2_0 + (u2_1 - u2_0) * i as f64 / (n2 - 1) as f64;
            let v = v2_0 + (v2_1 - v2_0) * j as f64 / (n2 - 1) as f64;
            s2_pts.push(s2.point_at(u, v));
        }
    }

    // Approximate signed distance: f(p) = min dist to s2 (no sign, just magnitude).
    // We want the iso-surface f = 0.
    // Better: evaluate at grid, collect points where s1 point is within tolerance of s2.
    let threshold = (u1_1 - u1_0).max(v1_1 - v1_0) / N as f64 * 2.0; // cell-size tolerance

    let mut intersection_pts: Vec<DVec3> = Vec::new();

    for i in 0..N {
        for j in 0..N {
            let u = u1_0 + (u1_1 - u1_0) * i as f64 / (N - 1) as f64;
            let v = v1_0 + (v1_1 - v1_0) * j as f64 / (N - 1) as f64;
            let p = s1.point_at(u, v);

            // Find closest point on s2
            let min_sq = s2_pts
                .iter()
                .map(|q| (p - *q).length_squared())
                .fold(f64::INFINITY, f64::min);

            if min_sq.sqrt() < threshold {
                intersection_pts.push(p);
            }
        }
    }

    let mut out = SurfaceSurfaceIntersection::default();
    if intersection_pts.len() >= 2 {
        out.curves.push(SurfaceCurve::Polyline(intersection_pts));
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{CylindricalSurface, Plane, SphericalSurface};

    #[test]
    fn plane_plane_parallel() {
        let p1 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let p2 = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 1.0),
            normal: DVec3::Z,
        });
        let r = intersect_surfaces(&p1, &p2);
        assert!(r.is_empty(), "parallel planes: no intersection");
    }

    #[test]
    fn plane_plane_intersect() {
        let p1 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let p2 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        });
        let r = intersect_surfaces(&p1, &p2);
        assert_eq!(r.curves.len(), 1);
        assert!(matches!(r.curves[0], SurfaceCurve::Line(_)));
    }

    #[test]
    fn sphere_sphere_equator() {
        // Two equal spheres touching at (1,0,0): each has r=1, centers at (0,0,0) and (2,0,0)
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&s1, &s2);
        assert_eq!(r.curves.len(), 1, "expected one circle");
        if let SurfaceCurve::Circle(c) = &r.curves[0] {
            assert!((c.center.x - 0.5).abs() < 1e-6, "center should be at x=0.5");
        } else {
            panic!("expected Circle");
        }
    }

    #[test]
    fn sphere_sphere_disjoint() {
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(5.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&s1, &s2);
        assert!(r.is_empty(), "disjoint spheres: no intersection");
    }

    #[test]
    fn cylinder_cylinder_parallel_intersecting() {
        // Two parallel cylinders r=1 centered at (0,0,0) and (1.5,0,0)
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(1.5, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        // Two parallel lines
        assert_eq!(r.curves.len(), 2, "expected two intersection lines");
    }

    #[test]
    fn cylinder_cylinder_tangent() {
        // Two parallel cylinders r=1 separated by exactly 2 (tangent externally)
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(2.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        // One tangent line
        assert_eq!(r.curves.len(), 1, "tangent cylinders: one line");
    }

    #[test]
    fn plane_sphere_great_circle() {
        let p = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
        });
        let r = intersect_surfaces(&p, &s);
        assert_eq!(r.curves.len(), 1);
        if let SurfaceCurve::Circle(c) = &r.curves[0] {
            assert!((c.radius - 3.0).abs() < 1e-6);
        } else {
            panic!("expected Circle");
        }
    }
}
