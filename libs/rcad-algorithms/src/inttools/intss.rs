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

use glam::DVec3;
use rcad_kernel::geom::{
    Circle3, ConicalSurface, Curve2d, Curve3, CylindricalSurface, Ellipse3, Hyperbola3, Line3,
    Parabola3, Plane, SphericalSurface, Surface3, SurfaceEval, any_perpendicular,
};

use crate::inttools::{
    pcurve_derive::{
        circle_pcurve_on_cylinder, circle_pcurve_on_plane, circle_pcurve_on_sphere,
        ellipse_pcurve_on_plane, fallback_pcurve_by_projection, line_pcurve_on_cylinder,
        line_pcurve_on_plane, polyline_pcurve_by_projection,
    },
    plane_cone::{PlaneConicalResult, intersect_plane_cone},
    plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder},
    plane_plane::{PlanePlaneResult, intersect_plane_plane},
    plane_sphere::{PlaneSphereResult, intersect_plane_sphere},
};
use crate::tolerance::{TOLERANCE_ABS, vectors_parallel};

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
    /// An exact analytic parabola.
    Parabola(Parabola3),
    /// An exact analytic hyperbola (single branch representation).
    Hyperbola(Hyperbola3),
    /// A tangent point (zero-dimensional contact).
    Point(DVec3),
    /// Numerically sampled polyline (fallback for non-analytic pairs).
    Polyline(Vec<DVec3>),
}

/// One intersection result: 3D curve plus optional PCurves on each surface.
#[derive(Debug, Clone)]
pub struct SurfaceIntersectionResult {
    pub curve_3d: SurfaceCurve,
    /// PCurve on surface A (populated in Task 3+).
    pub pcurve_on_a: Option<Curve2d>,
    /// PCurve on surface B (populated in Task 3+).
    pub pcurve_on_b: Option<Curve2d>,
}

/// All intersection curves / components found between two surfaces.
#[derive(Debug, Clone, Default)]
pub struct SurfaceSurfaceIntersection {
    pub curves: Vec<SurfaceIntersectionResult>,
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
    intersect_surfaces_with_density(s1, s2, 48)
}

/// Like [`intersect_surfaces`] but lets the caller specify the grid density `n`
/// for the numerical fallback.  Analytic pairs always return exact results
/// regardless of `n`.  The numerical fallback uses an `n×n` parameter-space
/// grid to find sign-change crossings.
///
/// Higher `n` gives more accurate intersection polylines at the cost of O(n²)
/// work.  The default used by [`intersect_surfaces`] is `n = 48`.
pub fn intersect_surfaces_with_density(
    s1: &Surface3,
    s2: &Surface3,
    grid_n: usize,
) -> SurfaceSurfaceIntersection {
    use Surface3::*;
    match (s1, s2) {
        // ── Plane × * ─────────────────────────────────────────────────────
        (Plane(p1), Plane(p2)) => plane_x_plane(p1, p2),
        (Plane(p), Sphere(s)) | (Sphere(s), Plane(p)) => plane_x_sphere(p, s),
        (Plane(p), Cylinder(c)) | (Cylinder(c), Plane(p)) => plane_x_cylinder(p, c),
        (Plane(p), Cone(c)) | (Cone(c), Plane(p)) => plane_x_cone(p, c),
        (Plane(p), Torus(t)) | (Torus(t), Plane(p)) => torus_x_plane(t, p),

        // ── Sphere × * ────────────────────────────────────────────────────
        (Sphere(s1), Sphere(s2)) => sphere_x_sphere(s1, s2),
        (Sphere(s), Cylinder(c)) | (Cylinder(c), Sphere(s)) => sphere_x_cylinder(s, c),
        (Sphere(s), Cone(c)) | (Cone(c), Sphere(s)) => sphere_x_cone(s, c),

        // ── Cylinder × Cylinder ───────────────────────────────────────────
        (Cylinder(c1), Cylinder(c2)) => cylinder_x_cylinder(c1, c2),

        // ── Cylinder × Cone and Cone × Cone fall through to numeric ───────
        // ── All others → numeric marching ─────────────────────────────────
        _ => numeric_intss_with_density(s1, s2, grid_n),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Plane × Plane
// ──────────────────────────────────────────────────────────────────────────────

fn plane_x_plane(p1: &Plane, p2: &Plane) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_plane(p1, p2) {
        PlanePlaneResult::Line(l) => {
            let pca = line_pcurve_on_plane(&l, p1);
            let pcb = line_pcurve_on_plane(&l, p2);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
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
            let pca = circle_pcurve_on_plane(&c, p);
            let pcb = circle_pcurve_on_sphere(&c, s);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(c),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneSphereResult::TangentPoint(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        PlaneSphereResult::NoIntersection => {}
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Plane × Cylinder
// ──────────────────────────────────────────────────────────────────────────────

fn plane_x_cylinder(p: &Plane, c: &CylindricalSurface) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_cylinder(p, c) {
        PlaneCylinderResult::Circle(circ) => {
            let pca = circle_pcurve_on_plane(&circ, p);
            let pcb = circle_pcurve_on_cylinder(&circ, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneCylinderResult::Ellipse(e) => {
            let pca = ellipse_pcurve_on_plane(&e, p);
            let pcb = fallback_pcurve_by_projection(
                &Curve3::Ellipse(e),
                &[0.0, TAU],
                &Surface3::Cylinder(*c),
            );
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Ellipse(e),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneCylinderResult::TangentLine(l) => {
            let pca = line_pcurve_on_plane(&l, p);
            let pcb = line_pcurve_on_cylinder(&l, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneCylinderResult::TwoLines(l1, l2) => {
            let pca1 = line_pcurve_on_plane(&l1, p);
            let pcb1 = line_pcurve_on_cylinder(&l1, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l1),
                pcurve_on_a: Some(pca1),
                pcurve_on_b: Some(pcb1),
            });
            let pca2 = line_pcurve_on_plane(&l2, p);
            let pcb2 = line_pcurve_on_cylinder(&l2, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l2),
                pcurve_on_a: Some(pca2),
                pcurve_on_b: Some(pcb2),
            });
        }
        PlaneCylinderResult::NoIntersection => {}
    }
    out
}

// ──────────────────────────────────────────────────────────────────────────────
// Plane × Cone
// ──────────────────────────────────────────────────────────────────────────────

fn plane_x_cone(p: &Plane, c: &ConicalSurface) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    let cone_surf = Surface3::Cone(*c);
    match intersect_plane_cone(p, c) {
        PlaneConicalResult::Circle(circ) => {
            let pca = circle_pcurve_on_plane(&circ, p);
            let pcb = fallback_pcurve_by_projection(&Curve3::Circle(circ), &[0.0, TAU], &cone_surf);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::Ellipse(e) => {
            let pca = ellipse_pcurve_on_plane(&e, p);
            let pcb = fallback_pcurve_by_projection(&Curve3::Ellipse(e), &[0.0, TAU], &cone_surf);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Ellipse(e),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::SingleLine(l) => {
            let pca = line_pcurve_on_plane(&l, p);
            let pcb = fallback_pcurve_by_projection(&Curve3::Line(l), &[-10.0, 10.0], &cone_surf);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::TwoLines(l1, l2) => {
            let pca1 = line_pcurve_on_plane(&l1, p);
            let pcb1 = fallback_pcurve_by_projection(&Curve3::Line(l1), &[-10.0, 10.0], &cone_surf);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l1),
                pcurve_on_a: Some(pca1),
                pcurve_on_b: Some(pcb1),
            });
            let pca2 = line_pcurve_on_plane(&l2, p);
            let pcb2 = fallback_pcurve_by_projection(&Curve3::Line(l2), &[-10.0, 10.0], &cone_surf);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l2),
                pcurve_on_a: Some(pca2),
                pcurve_on_b: Some(pcb2),
            });
        }
        PlaneConicalResult::Parabola(par) => {
            // Sample over a reasonable bounded domain for the PCurves
            let pca = fallback_pcurve_by_projection(&Curve3::Parabola(par), &[-20.0, 20.0], &Surface3::Plane(*p));
            let pcb = fallback_pcurve_by_projection(&Curve3::Parabola(par), &[-20.0, 20.0], &cone_surf);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Parabola(par),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::Hyperbola(hyp) => {
            // Each branch sampled separately; use the principal branch domain
            let pca = fallback_pcurve_by_projection(&Curve3::Hyperbola(hyp), &[-10.0, 10.0], &Surface3::Plane(*p));
            let pcb = fallback_pcurve_by_projection(&Curve3::Hyperbola(hyp), &[-10.0, 10.0], &cone_surf);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Hyperbola(hyp),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::Point(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
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
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Point(center),
            pcurve_on_a: None,
            pcurve_on_b: None,
        });
    } else {
        let circle = Circle3 {
            center,
            normal: axis,
            radius: r_circle,
        };
        let pca = circle_pcurve_on_sphere(&circle, s1);
        let pcb = circle_pcurve_on_sphere(&circle, s2);
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Circle(circle),
            pcurve_on_a: Some(pca),
            pcurve_on_b: Some(pcb),
        });
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
            let circle = Circle3 {
                center: s.center,
                normal: c.axis,
                radius: c.radius,
            };
            let pca = circle_pcurve_on_sphere(&circle, s);
            let pcb = circle_pcurve_on_cylinder(&circle, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circle),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        } else {
            let dz = dz_sq.sqrt();
            for &sign in &[1.0f64, -1.0] {
                let center = s.center + c.axis * (sign * dz);
                let circle = Circle3 {
                    center,
                    normal: c.axis,
                    radius: c.radius,
                };
                let pca = circle_pcurve_on_sphere(&circle, s);
                let pcb = circle_pcurve_on_cylinder(&circle, c);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circle),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
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
                    let circle = Circle3 {
                        center,
                        normal: c.axis,
                        radius: r_sol,
                    };
                    let pca = circle_pcurve_on_sphere(&circle, s);
                    let pcb = fallback_pcurve_by_projection(
                        &Curve3::Circle(circle),
                        &[0.0, std::f64::consts::TAU],
                        &Surface3::Cone(*c),
                    );
                    out.curves.push(SurfaceIntersectionResult {
                        curve_3d: SurfaceCurve::Circle(circle),
                        pcurve_on_a: Some(pca),
                        pcurve_on_b: Some(pcb),
                    });
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
            let line = Line3 {
                origin: pt,
                direction: c1.axis,
            };
            let pca = line_pcurve_on_cylinder(&line, c1);
            let pcb = line_pcurve_on_cylinder(&line, c2);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(line),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
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
// Torus × Plane
// ──────────────────────────────────────────────────────────────────────────────

/// Torus-plane intersection.
///
/// **Analytic case — plane ⊥ torus axis**:
///   If `plane.normal ∥ torus.axis` (the plane is perpendicular to the torus
///   axis), the intersection consists of up to two circles coaxial with the
///   torus.  Let `d` be the signed distance from the torus center to the
///   plane along the axis.  The plane cuts the torus tube at the radii
///   `sqrt((R ± sqrt(r²-d²))²)` where R is the major radius and r the minor.
///   More simply: the intersection circles have radii `R + sqrt(r²-d²)` and
///   `R - sqrt(r²-d²)` (the latter exists only when it is positive).
///   This simplifies to: for d² ≤ r², the two circle radii are
///   `R ± sqrt(r² - d²)`, both centered at the torus center projected onto
///   the plane.
///
/// **All other planes** fall back to numerical marching.
fn torus_x_plane(
    torus: &rcad_kernel::geom::ToroidalSurface,
    plane: &rcad_kernel::geom::Plane,
) -> SurfaceSurfaceIntersection {
    let axis = torus.axis.normalize();
    let normal = plane.normal.normalize();

    // Check whether the plane is perpendicular to the torus axis.
    // cos(angle between normal and axis): perpendicular ⟺ |cos| ≈ 1.
    let cos_angle = axis.dot(normal).abs();
    const PERP_TOL: f64 = 1e-6;

    if (cos_angle - 1.0).abs() > PERP_TOL {
        // Not perpendicular — fall back to numerical
        return numeric_intss(
            &Surface3::Torus(*torus),
            &Surface3::Plane(*plane),
        );
    }

    // Signed distance from torus center to the plane along the axis.
    let d = (plane.origin - torus.center).dot(normal) * normal.dot(axis).signum();
    let d_sq = d * d;
    let r_sq = torus.minor_radius * torus.minor_radius;

    if d_sq > r_sq + TOLERANCE_ABS {
        // Plane misses the torus tube
        return SurfaceSurfaceIntersection::default();
    }

    // Two intersection circles (when d²=r² they merge into one)
    let delta = (r_sq - d_sq).max(0.0).sqrt();
    let r1 = torus.major_radius + delta;
    let r2 = torus.major_radius - delta;

    // Center of intersection circles: projection of torus center onto plane.
    let center_proj = torus.center + axis * d;

    // Build circle normal (same as plane normal, oriented to match plane)
    let circle_normal = normal;

    let mut out = SurfaceSurfaceIntersection::default();

    // Outer circle (r1 > 0 always for valid torus)
    if r1 > TOLERANCE_ABS {
        let pcurve_a = pcurve_for_torus_circle(torus, center_proj, r1, plane);
        let pcurve_b = crate::inttools::pcurve_derive::circle_pcurve_on_plane(
            &rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r1,
            },
            plane,
        );
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Circle(rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r1,
            }),
            pcurve_on_a: Some(pcurve_a),
            pcurve_on_b: Some(pcurve_b),
        });
    }

    // Inner circle (r2 > 0 only if delta < major_radius)
    if r2 > TOLERANCE_ABS {
        let pcurve_a = pcurve_for_torus_circle(torus, center_proj, r2, plane);
        let pcurve_b = crate::inttools::pcurve_derive::circle_pcurve_on_plane(
            &rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r2,
            },
            plane,
        );
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Circle(rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r2,
            }),
            pcurve_on_a: Some(pcurve_a),
            pcurve_on_b: Some(pcurve_b),
        });
    }

    out
}

/// Compute a UV PCurve for a circle of intersection on a torus.
///
/// The circle has given center and radius in 3D, lying on the plane
/// perpendicular to the torus axis. Each point on the circle corresponds to
/// a unique major angle u (azimuth around the torus axis) and a fixed minor
/// angle v (around the tube). Returns a BSpline2 approximation.
fn pcurve_for_torus_circle(
    torus: &rcad_kernel::geom::ToroidalSurface,
    circle_center: DVec3,
    circle_radius: f64,
    _plane: &rcad_kernel::geom::Plane,
) -> rcad_kernel::geom::Curve2d {
    use rcad_kernel::projection::closest_point_on_surface;

    // Sample the circle in 3D and project each point onto the torus UV domain.
    let n = 33_usize;
    let u_ax = rcad_kernel::geom::any_perpendicular(torus.axis);
    let v_ax = torus.axis.cross(u_ax).normalize();

    let mut uv_pts: Vec<glam::DVec2> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / (n - 1) as f64;
            let p = circle_center + (u_ax * theta.cos() + v_ax * theta.sin()) * circle_radius;
            let proj = closest_point_on_surface(&Surface3::Torus(*torus), p, 16);
            glam::DVec2::new(proj.params.0, proj.params.1)
        })
        .collect();

    // Unwrap u discontinuities across the 2π seam.
    for i in 1..uv_pts.len() {
        let du = uv_pts[i].x - uv_pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut uv_pts[i..] { p.x -= std::f64::consts::TAU; }
        } else if du < -std::f64::consts::PI {
            for p in &mut uv_pts[i..] { p.x += std::f64::consts::TAU; }
        }
    }

    rcad_kernel::fit::interpolate_points_2d(&uv_pts)
        .map(rcad_kernel::geom::Curve2d::BSpline)
        .unwrap_or_else(|_| rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: uv_pts[0],
            direction: glam::DVec2::X,
        }))
}



/// Numerical surface-surface intersection via sign-change edge marching.
///
/// **Algorithm**:
/// 1. Sample `s1` on an N×N grid; for each sample, compute approximate distance to `s2`
///    using a pre-sampled `s2` grid.
/// 2. Detect edges (horizontal or vertical between adjacent grid cells) where the
///    distance changes sign (one end < threshold, other ≥ threshold or vice versa).
///    Linearly interpolate each crossing → candidate intersection points.
/// 3. BFS-greedy sort: start from any unvisited point, repeatedly extend the chain
///    by picking the nearest unvisited neighbor. Repeat until all points are visited.
///    This produces ordered polylines suitable for UV splitting.
fn numeric_intss(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection {
    numeric_intss_with_density(s1, s2, 48)
}

/// Same as `numeric_intss` but with configurable grid density N.
pub fn numeric_intss_with_density(
    s1: &Surface3,
    s2: &Surface3,
    n: usize,
) -> SurfaceSurfaceIntersection {
    numeric_intss_impl(s1, s2, n, None, None)
}

/// Same as `numeric_intss_with_density` but uses caller-supplied UV domains
/// for s1 and s2 instead of `default_domain()`.  Pass `None` for either to
/// use the surface's own default domain (with infinite-domain clamping).
pub fn numeric_intss_with_domains(
    s1: &Surface3,
    s2: &Surface3,
    n: usize,
    dom1_override: Option<[f64; 4]>,
    dom2_override: Option<[f64; 4]>,
) -> SurfaceSurfaceIntersection {
    numeric_intss_impl(s1, s2, n, dom1_override, dom2_override)
}

fn numeric_intss_impl(
    s1: &Surface3,
    s2: &Surface3,
    n: usize,
    dom1_override: Option<[f64; 4]>,
    dom2_override: Option<[f64; 4]>,
) -> SurfaceSurfaceIntersection {
    let dom1 = s1.default_domain();
    let dom2 = s2.default_domain();

    // Clamp infinite domains. For cylinders the v-domain is [-∞, +∞]; we use
    // a range large enough to cover any practical intersection geometry.
    const DOMAIN_CLAMP: f64 = 10.0;
    let clamp_dom = |[u0, u1, v0, v1]: [f64; 4]| -> [f64; 4] {
        [
            if u0.is_finite() { u0 } else { -DOMAIN_CLAMP },
            if u1.is_finite() { u1 } else { DOMAIN_CLAMP },
            if v0.is_finite() { v0 } else { -DOMAIN_CLAMP },
            if v1.is_finite() { v1 } else { DOMAIN_CLAMP },
        ]
    };
    let [u1_0, u1_1, v1_0, v1_1] = dom1_override.unwrap_or_else(|| clamp_dom(dom1));
    let [u2_0, u2_1, v2_0, v2_1] = dom2_override.unwrap_or_else(|| clamp_dom(dom2));

    // Pre-sample s2 on a grid for fast approximate distance computation
    let n2 = n.min(48);
    let mut s2_pts: Vec<DVec3> = Vec::with_capacity(n2 * n2);
    for i in 0..n2 {
        for j in 0..n2 {
            let u = u2_0 + (u2_1 - u2_0) * i as f64 / (n2 - 1).max(1) as f64;
            let v = v2_0 + (v2_1 - v2_0) * j as f64 / (n2 - 1).max(1) as f64;
            s2_pts.push(s2.point_at(u, v));
        }
    }

    // Approximate distance from 3D point to s2 surface via closest sample
    let approx_dist_to_s2 = |p: DVec3| -> f64 {
        s2_pts
            .iter()
            .map(|q| (*q - p).length())
            .fold(f64::INFINITY, f64::min)
    };

    // Threshold: treated as "on the surface" if distance < this.
    // Use the average cell size on s1 as a reference scale.
    let du = (u1_1 - u1_0) / n as f64;
    let dv = (v1_1 - v1_0) / n as f64;
    let p00 = s1.point_at(u1_0, v1_0);
    let p10 = s1.point_at(u1_0 + du, v1_0);
    let p01 = s1.point_at(u1_0, v1_0 + dv);
    let cell_size = (p10 - p00).length().max((p01 - p00).length()).max(1e-6);
    let threshold = cell_size * 2.0;

    // Compute distance at each grid point
    let nn = n + 1; // grid has (n+1) × (n+1) nodes
    let mut dist: Vec<f64> = Vec::with_capacity(nn * nn);
    let mut pts: Vec<DVec3> = Vec::with_capacity(nn * nn);
    for i in 0..nn {
        for j in 0..nn {
            let u = u1_0 + (u1_1 - u1_0) * i as f64 / n as f64;
            let v = v1_0 + (v1_1 - v1_0) * j as f64 / n as f64;
            let p = s1.point_at(u, v);
            pts.push(p);
            dist.push(approx_dist_to_s2(p));
        }
    }

    let idx = |i: usize, j: usize| i * nn + j;

    // Find sign-change edges and interpolate crossing points
    let mut crossing_pts: Vec<DVec3> = Vec::new();

    // Horizontal edges: (i,j) — (i, j+1)
    for i in 0..nn {
        for j in 0..n {
            let a = idx(i, j);
            let b = idx(i, j + 1);
            let da = dist[a];
            let db = dist[b];
            // Sign change: one below threshold, other above (or vice versa)
            if (da < threshold) != (db < threshold) {
                let t = if (da - db).abs() < 1e-15 {
                    0.5
                } else {
                    (threshold - da) / (db - da)
                };
                let t = t.clamp(0.0, 1.0);
                crossing_pts.push(pts[a].lerp(pts[b], t));
            }
        }
    }

    // Vertical edges: (i,j) — (i+1, j)
    for i in 0..n {
        for j in 0..nn {
            let a = idx(i, j);
            let b = idx(i + 1, j);
            let da = dist[a];
            let db = dist[b];
            if (da < threshold) != (db < threshold) {
                let t = if (da - db).abs() < 1e-15 {
                    0.5
                } else {
                    (threshold - da) / (db - da)
                };
                let t = t.clamp(0.0, 1.0);
                crossing_pts.push(pts[a].lerp(pts[b], t));
            }
        }
    }

    let mut out = SurfaceSurfaceIntersection::default();

    if crossing_pts.len() < 2 {
        return out;
    }

    // BFS-greedy ordering: connect nearest unvisited neighbors into chains.
    // This works well for smooth curves; for self-intersecting surfaces it may
    // produce slightly wrong orderings near the crossing, which is acceptable
    // for topological boolean operations.
    let ordered = greedy_order_points(crossing_pts);

    for chain in ordered {
        if chain.len() < 2 {
            continue;
        }
        let pca = polyline_pcurve_by_projection(&chain, s1);
        let pcb = polyline_pcurve_by_projection(&chain, s2);
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Polyline(chain),
            pcurve_on_a: pca,
            pcurve_on_b: pcb,
        });
    }

    out
}

/// Greedy nearest-neighbor ordering of a point cloud into one or more chains.
///
/// Returns a `Vec` of chains (each chain is `Vec<DVec3>`).
/// Points that can't be extended within `gap_tol` start a new chain.
/// After chain formation, chains are stitched together when their endpoints
/// are close enough, producing fewer, longer chains (typically one closed loop
/// for a single intersection curve).
fn greedy_order_points(pts: Vec<DVec3>) -> Vec<Vec<DVec3>> {
    if pts.is_empty() {
        return vec![];
    }

    // Estimate gap tolerance from average nearest-neighbor distance
    // (rough: use 3x the median distance between sorted x-coordinates)
    let gap_tol = {
        let mut dists: Vec<f64> = Vec::with_capacity(pts.len());
        for i in 0..pts.len() {
            let mut best = f64::INFINITY;
            for j in 0..pts.len() {
                if i != j {
                    let d = (pts[i] - pts[j]).length();
                    if d < best {
                        best = d;
                    }
                }
            }
            dists.push(best);
        }
        dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = dists[dists.len() / 2];
        (median * 5.0).max(1e-9)
    };

    // Stitch gap tolerance: more generous than within-chain growth.
    // Allow up to 3× the within-chain gap to merge separate chains.
    let stitch_tol = gap_tol * 3.0;

    let mut used = vec![false; pts.len()];
    let mut chains: Vec<Vec<DVec3>> = Vec::new();

    loop {
        // Find first unused point
        let start = match used.iter().position(|&u| !u) {
            Some(i) => i,
            None => break,
        };
        used[start] = true;
        let mut chain = vec![pts[start]];

        loop {
            let last = *chain.last().expect("chain is non-empty (starts with 1 element)");
            // Find nearest unused point within gap_tol
            let mut best_dist = gap_tol;
            let mut best_idx = None;
            for (i, &used_i) in used.iter().enumerate() {
                if !used_i {
                    let d = (pts[i] - last).length();
                    if d < best_dist {
                        best_dist = d;
                        best_idx = Some(i);
                    }
                }
            }
            match best_idx {
                Some(idx) => {
                    used[idx] = true;
                    chain.push(pts[idx]);
                }
                None => break,
            }
        }

        chains.push(chain);
    }

    // ── Chain stitching ──────────────────────────────────────────────────────
    // Repeatedly merge pairs of chains whose endpoints are within stitch_tol.
    // This turns fragmented arc segments into a single closed loop.
    let mut changed = true;
    while changed && chains.len() > 1 {
        changed = false;
        'outer: for i in 0..chains.len() {
            for j in (i + 1)..chains.len() {
                let end_i = *chains[i].last().expect("chains[i] is non-empty");
                let start_j = chains[j][0];
                let end_j = *chains[j].last().expect("chains[j] is non-empty");
                let start_i = chains[i][0];

                // Determine merge direction
                let (merge_rev_j, close_enough) =
                    if (end_i - start_j).length() <= stitch_tol {
                        (false, true) // i + j
                    } else if (end_i - end_j).length() <= stitch_tol {
                        (true, true)  // i + reversed j
                    } else if (end_j - start_i).length() <= stitch_tol {
                        // j + i: handled next iteration via swapped roles
                        (false, false)
                    } else if (start_j - start_i).length() <= stitch_tol {
                        (false, false)
                    } else {
                        (false, false)
                    };

                if close_enough {
                    let chain_j = chains.remove(j);
                    let appended: Vec<DVec3> = if merge_rev_j {
                        chain_j.into_iter().rev().collect()
                    } else {
                        chain_j
                    };
                    chains[i].extend(appended);
                    changed = true;
                    break 'outer;
                }

                // Also handle: j ends near i start → prepend j to i
                if (end_j - start_i).length() <= stitch_tol {
                    let chain_j = chains.remove(j);
                    let mut merged = chain_j;
                    merged.extend(chains[i].drain(..));
                    chains[i] = merged;
                    changed = true;
                    break 'outer;
                }
                // j start near i start → prepend reversed j
                if (start_j - start_i).length() <= stitch_tol {
                    let chain_j = chains.remove(j);
                    let mut merged: Vec<DVec3> = chain_j.into_iter().rev().collect();
                    merged.extend(chains[i].drain(..));
                    chains[i] = merged;
                    changed = true;
                    break 'outer;
                }
            }
        }
    }

    chains
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
        assert!(matches!(r.curves[0].curve_3d, SurfaceCurve::Line(_)));
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
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
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
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            assert!((c.radius - 3.0).abs() < 1e-6);
        } else {
            panic!("expected Circle");
        }
    }

    #[test]
    fn cylinder_cylinder_perpendicular_steinmetz() {
        // Two perpendicular cylinders r=1 — Steinmetz configuration
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(0.0, -2.0, 0.0),
            axis: DVec3::Y,
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(-2.0, 0.0, 0.0),
            axis: DVec3::X,
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        // Should find the Steinmetz intersection curve(s)
        assert!(!r.curves.is_empty(), "expected at least one intersection curve, got none");
        // The Steinmetz intersection is one or two closed space curves
        if let SurfaceCurve::Polyline(pts) = &r.curves[0].curve_3d {
            assert!(pts.len() >= 4, "polyline should have ≥4 points, got {}", pts.len());
        }
    }

    #[test]
    fn torus_perpendicular_plane_gives_circles() {
        use rcad_kernel::geom::ToroidalSurface;

        // Torus with axis=Z, centered at origin, R=5, r=1.
        // Plane at z=0 (perpendicular to the axis) intersects the torus
        // in two concentric circles with radii R+r=6 and R-r=4.
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        });
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let r = intersect_surfaces(&torus, &plane);
        assert_eq!(
            r.curves.len(),
            2,
            "torus ∩ perp-plane should give 2 circles, got {}",
            r.curves.len()
        );

        // Collect radii
        let mut radii: Vec<f64> = r
            .curves
            .iter()
            .filter_map(|c| {
                if let SurfaceCurve::Circle(circ) = &c.curve_3d {
                    Some(circ.radius)
                } else {
                    None
                }
            })
            .collect();
        radii.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        assert_eq!(radii.len(), 2, "expected 2 Circle3 results");
        assert!(
            (radii[0] - 4.0).abs() < 1e-6,
            "inner circle radius should be 4, got {}",
            radii[0]
        );
        assert!(
            (radii[1] - 6.0).abs() < 1e-6,
            "outer circle radius should be 6, got {}",
            radii[1]
        );
    }
}

