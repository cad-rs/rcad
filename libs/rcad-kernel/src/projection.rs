//! Closest-point projection from a 3D point onto a curve or surface.
//!
//! Analogous to OCCT `GeomAPI_ProjectPointOnCurve` and
//! `GeomAPI_ProjectPointOnSurf`.
//!
//! # Strategy
//! - **Analytic surfaces** (Plane, Cylinder, Sphere, Cone, Torus): closed-form
//!   projection — fast and exact.
//! - **All curves** and **parametric surfaces** (BSpline, Bezier, Offset,
//!   LinearExtrusion, Revolution): sample the domain uniformly to find the best
//!   initial guess, then refine with Newton-Raphson minimisation of `|P(t) - Q|²`.

use glam::DVec3;

use crate::geom::{Curve3, CurveEval, Surface3, SurfaceEval};

// ─────────────────────────────────────────────────────────────────────────────
// Result types
// ─────────────────────────────────────────────────────────────────────────────

/// Result of projecting a point onto a curve.
#[derive(Debug, Clone)]
pub struct CurveProjection {
    /// Nearest point on the curve.
    pub point: DVec3,
    /// Curve parameter at the nearest point.
    pub param: f64,
    /// Distance from the query point to the curve.
    pub distance: f64,
}

/// Result of projecting a point onto a surface.
#[derive(Debug, Clone)]
pub struct SurfaceProjection {
    /// Nearest point on the surface.
    pub point: DVec3,
    /// Surface parameter (u, v) at the nearest point.
    pub params: (f64, f64),
    /// Distance from the query point to the surface.
    pub distance: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Curve projection
// ─────────────────────────────────────────────────────────────────────────────

/// Project the point `query` onto `curve`, returning the nearest point on the
/// curve, its parameter value, and the Euclidean distance.
///
/// The curve is evaluated over its natural domain ([`CurveEval::default_domain`]).
///
/// # Algorithm
/// 1. Sample `n_samples` points uniformly in the domain.
/// 2. Take the sample with the smallest distance as the initial guess.
/// 3. Refine with Newton iterations minimising `f(t) = |C(t) - Q|²`:
///    `t_{i+1} = t_i - (C(t) - Q) · T(t) / (|T|² + (C - Q) · T')` where
///    `T = C'(t)` (approximated by finite difference).
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Circle3};
/// use rcad_kernel::projection::closest_point_on_curve;
///
/// let circle = Curve3::Circle(Circle3 {
///     center: DVec3::ZERO,
///     normal: DVec3::Z,
///     radius: 1.0,
/// });
/// let q = DVec3::new(2.0, 0.0, 0.0);
/// let result = closest_point_on_curve(&circle, q, 64);
/// assert!((result.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
/// ```
pub fn closest_point_on_curve(curve: &Curve3, query: DVec3, n_samples: usize) -> CurveProjection {
    let [t0_raw, t1_raw] = curve.default_domain();
    let n = n_samples.max(4);

    // For infinite domains (lines), use a heuristic finite sampling range
    // centered on the closest parameter analytically (dot product for lines).
    let (t0, t1) = if t0_raw.is_infinite() || t1_raw.is_infinite() {
        // Use the analytical projection for the domain center estimate
        let t_center = match curve {
            Curve3::Line(l) => (query - l.origin).dot(l.direction),
            _ => 0.0,
        };
        let span = 100.0_f64; // generous range around t_center
        (t_center - span, t_center + span)
    } else {
        (t0_raw, t1_raw)
    };

    // Step 1: coarse sampling
    let (mut best_t, mut best_dist) = {
        let mut bd = f64::INFINITY;
        let mut bt = t0;
        for i in 0..=n {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            let p = curve.point_at(t);
            let d = (p - query).length();
            if d < bd {
                bd = d;
                bt = t;
            }
        }
        (bt, bd)
    };

    // Step 2: Newton refinement
    // For infinite domains, don't clamp the Newton step
    let clamp_t = |t: f64| {
        if t0_raw.is_infinite() || t1_raw.is_infinite() {
            t
        } else {
            t.clamp(t0, t1)
        }
    };
    let dt = if (t1 - t0).is_finite() {
        (t1 - t0) * 1e-6
    } else {
        1e-6
    };
    for _ in 0..30 {
        let p = curve.point_at(best_t);
        let diff = p - query;
        // Finite-difference tangent
        let t_plus = best_t + dt;
        let t_minus = best_t - dt;
        let span = t_plus - t_minus;
        if span.abs() < 1e-20 {
            break;
        }
        let tangent = (curve.point_at(t_plus) - curve.point_at(t_minus)) / span;
        let tang_sq = tangent.dot(tangent);
        if tang_sq < 1e-20 {
            break;
        }
        // Second-order term (curvature denominator term)
        let curvature_approx = (curve.point_at(best_t + 2.0 * dt) - 2.0 * p
            + curve.point_at(best_t - 2.0 * dt))
            / (dt * dt);
        let denom = tang_sq + diff.dot(curvature_approx);
        let delta = diff.dot(tangent) / if denom.abs() > 1e-20 { denom } else { tang_sq };
        let new_t = clamp_t(best_t - delta);
        let new_dist = (curve.point_at(new_t) - query).length();
        if new_dist < best_dist {
            best_dist = new_dist;
            best_t = new_t;
        }
        if delta.abs() < 1e-10 {
            break;
        }
    }

    let best_point = curve.point_at(best_t);
    CurveProjection {
        point: best_point,
        param: best_t,
        distance: (best_point - query).length(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface projection
// ─────────────────────────────────────────────────────────────────────────────

/// Project the point `query` onto `surface`, returning the nearest point on the
/// surface, its (u, v) parameters, and the Euclidean distance.
///
/// Analytic surfaces (Plane, Sphere, Cylinder, Cone, Torus) use closed-form
/// formulae.  All other surfaces fall back to numerical sampling + Newton
/// refinement.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Surface3, SphericalSurface};
/// use rcad_kernel::projection::closest_point_on_surface;
///
/// let sphere = Surface3::Sphere(SphericalSurface {
///     center: DVec3::ZERO,
///     axis: DVec3::Z,
///     radius: 1.0,
/// });
/// let q = DVec3::new(3.0, 0.0, 0.0);
/// let result = closest_point_on_surface(&sphere, q, 16);
/// assert!((result.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
/// ```
pub fn closest_point_on_surface(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    use crate::geom::*;

    match surface {
        // ── Analytic closed-form projections ──────────────────────────────────
        Surface3::Plane(plane) => {
            // Project onto the infinite plane, no clamping needed.
            let d = (query - plane.origin).dot(plane.normal);
            let point = query - plane.normal * d;
            let distance = d.abs();
            // (u, v) = coordinates in plane (not strictly needed, but provide them)
            let u_axis = any_perpendicular(plane.normal);
            let v_axis = plane.normal.cross(u_axis);
            let diff = point - plane.origin;
            SurfaceProjection {
                point,
                params: (diff.dot(u_axis), diff.dot(v_axis)),
                distance,
            }
        }

        Surface3::Sphere(sph) => {
            let v = query - sph.center;
            let len = v.length();
            let point = if len < 1e-14 {
                sph.center + sph.radius * DVec3::X // degenerate: pick arbitrary
            } else {
                sph.center + v / len * sph.radius
            };
            // Compute (theta, phi) from point relative to center
            let w = (point - sph.center).normalize_or_zero();
            let u_axis = any_perpendicular(sph.axis);
            let v_axis = sph.axis.cross(u_axis);
            let theta = w.dot(sph.axis).clamp(-1.0, 1.0).acos();
            let phi = w.dot(v_axis).atan2(w.dot(u_axis));
            SurfaceProjection {
                point,
                params: (phi, theta),
                distance: (point - query).length(),
            }
        }

        Surface3::Cylinder(cyl) => {
            // Project by collapsing along axis, then normalizing radial component.
            let v = query - cyl.origin;
            let along = v.dot(cyl.axis);
            let radial = v - cyl.axis * along;
            let radial_len = radial.length();
            let point = if radial_len < 1e-14 {
                cyl.origin + cyl.axis * along + cyl.radius * any_perpendicular(cyl.axis)
            } else {
                cyl.origin + cyl.axis * along + radial / radial_len * cyl.radius
            };
            let u_axis = any_perpendicular(cyl.axis);
            let v_axis = cyl.axis.cross(u_axis);
            let r = (point - cyl.origin - cyl.axis * along).normalize_or_zero();
            let theta = r.dot(v_axis).atan2(r.dot(u_axis));
            SurfaceProjection {
                point,
                params: (theta, along),
                distance: (point - query).length(),
            }
        }

        Surface3::Cone(cone) => {
            // Project onto the cone: find the closest generator line.
            let v = query - cone.apex;
            let along = v.dot(cone.axis);
            let radial = v - cone.axis * along;
            let radial_len = radial.length();
            let half = cone.half_angle_rad;
            // The foot on the cone satisfies r = s·tan(half), z = s
            // Minimize |Q - (apex + s*axis + s*tan(half)*r_hat)|²
            // → s = (along + radial_len*tan(half)) / (1 + tan(half)²)
            let tan_h = half.tan();
            let s = (along + radial_len * tan_h) / (1.0 + tan_h * tan_h);
            let s = s.max(0.0);
            let r_hat = if radial_len < 1e-14 {
                any_perpendicular(cone.axis)
            } else {
                radial / radial_len
            };
            let point = cone.apex + cone.axis * s + r_hat * s * tan_h;
            SurfaceProjection {
                point,
                params: (
                    s,
                    r_hat
                        .dot(any_perpendicular(cone.axis))
                        .atan2(r_hat.dot(cone.axis.cross(any_perpendicular(cone.axis)))),
                ),
                distance: (point - query).length(),
            }
        }

        Surface3::Torus(torus) => {
            // Step 1: project onto the major-radius circle in the equatorial plane.
            let v = query - torus.center;
            let along = v.dot(torus.axis);
            let radial = v - torus.axis * along;
            let radial_len = radial.length();
            let major_dir = if radial_len < 1e-14 {
                any_perpendicular(torus.axis)
            } else {
                radial / radial_len
            };
            let tube_center = torus.center + major_dir * torus.major_radius;
            // Step 2: project onto the tube circle.
            let w = query - tube_center;
            let w_len = w.length();
            let point = if w_len < 1e-14 {
                tube_center + major_dir * torus.minor_radius
            } else {
                tube_center + w / w_len * torus.minor_radius
            };
            let u = major_dir
                .dot(any_perpendicular(torus.axis))
                .atan2(major_dir.dot(torus.axis.cross(any_perpendicular(torus.axis))));
            let w_dir = (point - tube_center).normalize_or_zero();
            let v_param = w_dir.dot(torus.axis).atan2(w_dir.dot(major_dir));
            SurfaceProjection {
                point,
                params: (u, v_param),
                distance: (point - query).length(),
            }
        }

        // ── Numerical fallback for parametric surfaces ─────────────────────────
        _ => numeric_surface_projection(surface, query, n_samples),
    }
}

/// Numerical closest-point on a parametric surface via uniform sampling +
/// Newton refinement of `f(u,v) = |S(u,v) - Q|²`.
fn numeric_surface_projection(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    let [u0, u1, v0, v1] = surface.default_domain();
    let n = n_samples.max(4);

    // Coarse sampling
    let (mut best_u, mut best_v, mut best_dist) = {
        let mut bd = f64::INFINITY;
        let (mut bu, mut bv) = (u0, v0);
        for i in 0..=n {
            for j in 0..=n {
                let u = u0 + (u1 - u0) * i as f64 / n as f64;
                let v = v0 + (v1 - v0) * j as f64 / n as f64;
                let p = surface.point_at(u, v);
                let d = (p - query).length_squared();
                if d < bd {
                    bd = d;
                    bu = u;
                    bv = v;
                }
            }
        }
        (bu, bv, bd.sqrt())
    };

    // Newton refinement: gradient of ½|S(u,v)-Q|²
    let eps = ((u1 - u0) + (v1 - v0)) * 1e-6;
    for _ in 0..40 {
        let p = surface.point_at(best_u, best_v);
        let diff = p - query;
        // Partial derivatives via finite difference
        let pu = surface.point_at((best_u + eps).min(u1), best_v);
        let pum = surface.point_at((best_u - eps).max(u0), best_v);
        let pv = surface.point_at(best_u, (best_v + eps).min(v1));
        let pvm = surface.point_at(best_u, (best_v - eps).max(v0));
        let du = (pu - pum) / (2.0 * eps.min((best_u + eps).min(u1) - (best_u - eps).max(u0)));
        let dv = (pv - pvm) / (2.0 * eps.min((best_v + eps).min(v1) - (best_v - eps).max(v0)));
        // Gradient components: ∂f/∂u = diff · du, ∂f/∂v = diff · dv
        let gu = diff.dot(du);
        let gv = diff.dot(dv);
        // Hessian diagonal approximation (Gauss-Newton)
        let huu = du.dot(du);
        let hvv = dv.dot(dv);
        let huv = du.dot(dv);
        let det = huu * hvv - huv * huv;
        if det.abs() < 1e-20 {
            break;
        }
        let delta_u = (hvv * gu - huv * gv) / det;
        let delta_v = (huu * gv - huv * gu) / det;
        let new_u = (best_u - delta_u).clamp(u0, u1);
        let new_v = (best_v - delta_v).clamp(v0, v1);
        let new_dist = (surface.point_at(new_u, new_v) - query).length();
        if new_dist < best_dist {
            best_dist = new_dist;
            best_u = new_u;
            best_v = new_v;
        }
        if delta_u.abs() < 1e-10 && delta_v.abs() < 1e-10 {
            break;
        }
    }

    let best_point = surface.point_at(best_u, best_v);
    SurfaceProjection {
        point: best_point,
        params: (best_u, best_v),
        distance: (best_point - query).length(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn project_onto_plane() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        });
        let q = DVec3::new(3.0, 5.0, -2.0);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!(
            (r.point - DVec3::new(3.0, 0.0, -2.0)).length() < 1e-9,
            "point={}",
            r.point
        );
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        let q = DVec3::new(5.0, 0.0, 0.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        assert!((r.point - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 3.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_cylinder() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });
        let q = DVec3::new(3.0, 2.0, 0.0);
        let r = closest_point_on_surface(&cyl, q, 16);
        assert!((r.point - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 2.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_torus() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        // Point far along +X axis → nearest point on outer equator
        let q = DVec3::new(10.0, 0.0, 0.0);
        let r = closest_point_on_surface(&torus, q, 16);
        // Nearest should be at (4, 0, 0)
        assert!((r.point - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn project_onto_circle_curve() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let q = DVec3::new(2.0, 0.0, 0.0);
        let r = closest_point_on_curve(&circle, q, 64);
        assert!(
            (r.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6,
            "expected (1,0,0) got {}",
            r.point
        );
        assert!((r.distance - 1.0).abs() < 1e-6);
    }

    #[test]
    fn project_onto_line_curve() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        // Line has infinite domain; nearest point to (3, 4, 0) is (3, 0, 0)
        let q = DVec3::new(3.0, 4.0, 0.0);
        let r = closest_point_on_curve(&line, q, 32);
        let expected = DVec3::new(3.0, 0.0, 0.0);
        assert!(
            (r.point - expected).length() < 1e-4,
            "expected {:?} got {}",
            expected,
            r.point
        );
        assert!((r.distance - 4.0).abs() < 1e-4, "distance={}", r.distance);
    }

    #[test]
    fn project_onto_bspline_surface() {
        // Flat BSpline surface at z=0 over [0,1]²
        use crate::geom::BSplineSurface;
        let surf = Surface3::BSpline(BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0; 2]; 2],
        });
        let q = DVec3::new(0.5, 0.5, 5.0);
        let r = closest_point_on_surface(&surf, q, 8);
        assert!(
            (r.point - DVec3::new(0.5, 0.5, 0.0)).length() < 1e-4,
            "got {}",
            r.point
        );
    }
}
