//! Analytic intersection of an Ellipse3 with analytic surfaces.
//!
//! OCCT reference: IntAna_IntConicQuad only handles conic × Plane.
//! OCCT dispatches Ellipse × {Cylinder, Cone, Sphere} through the generic
//! numerical Edge-Face intersection (IntTools_EdgeFace), same as rcad's
//! ImplicitSurface numeric path.
//!
//! Ellipse × Plane: analytic A + B·cosθ + C·sinθ = 0 solve, matching OCCT.
//! Ellipse × {Cylinder, Cone, Sphere}: numeric Newton refinement, matching
//!    OCCT's generic EF fallback path in IntTools_EdgeFace.

use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// Hit from an ellipse-surface intersection.
pub struct EllipseSurfaceHit {
    pub point: DVec3,
    /// Parametric value on the ellipse (angle in radians).
    pub ellipse_param: f64,
}

/// Intersect an ellipse arc with a plane.
///
/// OCCT reference: IntAna_IntConicQuad — the ellipse is a conic (Geom_Conic)
/// evaluated at parameter θ. Substituting P(θ) = center + a·cosθ·u + b·sinθ·v
/// into the plane implicit n·(X - origin) = 0 gives:
///
///   n·(center - origin) + a·(n·u)·cosθ + b·(n·v)·sinθ = 0
///
/// i.e. A + B·cosθ + C·sinθ = 0 → R·cos(θ - φ) = -A
///
/// Returns 0, 1, or 2 hits within `t_range`.
///
/// same analytic method as IntAna_IntConicQuad for planes.
pub fn intersect_ellipse_plane(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    plane: &Plane,
) -> Vec<EllipseSurfaceHit> {
    intersect_ellipse_plane_with_tol(ellipse, t_range, plane, TOLERANCE_ABS)
}

/// Same as [`intersect_ellipse_plane`] with angular / amplitude margins from `geom_tol`.
pub fn intersect_ellipse_plane_with_tol(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    plane: &Plane,
    geom_tol: f64,
) -> Vec<EllipseSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);

    // Local orthonormal frame: u = major_dir, v = normal × major_dir
    let u = ellipse.major_dir.normalize();
    let v = ellipse.normal.cross(u).normalize();

    // Plane implicit: (P - plane.origin)·n = 0
    let d = (ellipse.center - plane.origin).dot(plane.normal);
    let a_coeff = ellipse.major_radius * u.dot(plane.normal);
    let b_coeff = ellipse.minor_radius * v.dot(plane.normal);

    // A + B·cos(θ) + C·sin(θ) = 0
    // → R·cos(θ - φ) = -A  where R = sqrt(B² + C²), φ = atan2(C, B)
    let r_amp = (a_coeff * a_coeff + b_coeff * b_coeff).sqrt();
    if r_amp < TOLERANCE_ABS {
        return vec![]; // ellipse parallel to plane
    }

    let ratio = -d / r_amp;
    if ratio.abs() > 1.0 + eps {
        return vec![];
    }
    let ratio = ratio.clamp(-1.0, 1.0);

    let phi = b_coeff.atan2(a_coeff);
    let alpha = ratio.acos();

    let [t0, t1] = t_range;
    let two_pi = 2.0 * std::f64::consts::PI;

    let mut hits = Vec::new();
    for theta in [phi + alpha, phi - alpha] {
        // Normalize theta to [0, 2π)
        let theta = ((theta % two_pi) + two_pi) % two_pi;

        if theta >= t0 - eps && theta <= t1 + eps {
            // Deduplicate: tangent case produces theta = phi + alpha = phi - alpha
            if hits.iter().any(|h: &EllipseSurfaceHit| (h.ellipse_param - theta).abs() < eps) {
                continue;
            }
            let point = ellipse.center
                + u * (ellipse.major_radius * theta.cos())
                + v * (ellipse.minor_radius * theta.sin());
            hits.push(EllipseSurfaceHit {
                point,
                ellipse_param: theta,
            });
        }
    }
    hits
}

/// Intersect an ellipse arc with a cylindrical surface.
///
/// OCCT IntAna_IntConicQuad has no ellipse×cylinder path — falls back to
/// IntTools_EdgeFace numeric.  rcad: Newton refinement on implicit equation.
/// OCCT-equivalent: numeric fallback matches OCCT's generic EF path.
pub fn intersect_ellipse_cylinder(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
) -> Vec<EllipseSurfaceHit> {
    intersect_ellipse_cylinder_with_tol(ellipse, t_range, cyl, TOLERANCE_ABS)
}

/// Same as [`intersect_ellipse_cylinder`] with Newton / range margins from `geom_tol`.
pub fn intersect_ellipse_cylinder_with_tol(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
    geom_tol: f64,
) -> Vec<EllipseSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    ellipse_vs_implicit_surface(
        ellipse,
        t_range,
        |p: DVec3| -> f64 {
            let v = p - cyl.origin;
            let along = v.dot(cyl.axis);
            let perp = v - cyl.axis * along;
            perp.length_squared() - cyl.radius * cyl.radius
        },
        eps,
    )
}

/// Intersect an ellipse arc with a spherical surface.
/// OCCT-equivalent: OCCT has no analytic path; uses generic EF numeric.
pub fn intersect_ellipse_sphere(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
) -> Vec<EllipseSurfaceHit> {
    intersect_ellipse_sphere_with_tol(ellipse, t_range, sph, TOLERANCE_ABS)
}

/// Same as [`intersect_ellipse_sphere`] with Newton / range margins from `geom_tol`.
pub fn intersect_ellipse_sphere_with_tol(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
    geom_tol: f64,
) -> Vec<EllipseSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    ellipse_vs_implicit_surface(
        ellipse,
        t_range,
        |p: DVec3| -> f64 { (p - sph.center).length_squared() - sph.radius * sph.radius },
        eps,
    )
}

/// Intersect an ellipse arc with a conical surface.
/// OCCT-equivalent: OCCT has no analytic path; uses generic EF numeric.
pub fn intersect_ellipse_cone(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
) -> Vec<EllipseSurfaceHit> {
    intersect_ellipse_cone_with_tol(ellipse, t_range, cone, TOLERANCE_ABS)
}

/// Same as [`intersect_ellipse_cone`] with Newton / range margins from `geom_tol`.
pub fn intersect_ellipse_cone_with_tol(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
    geom_tol: f64,
) -> Vec<EllipseSurfaceHit> {
    let cos2 = cone.half_angle_rad.cos().powi(2);
    let apex = cone.apex_point();
    let axis = cone.axis_dir();
    let eps = geom_tol.max(TOLERANCE_ABS);
    ellipse_vs_implicit_surface(
        ellipse,
        t_range,
        |p: DVec3| -> f64 {
            let v = p - apex;
            let along = v.dot(axis);
            let along2 = along * along;
            let len2 = v.length_squared();
            along2 - cos2 * len2
        },
        eps,
    )
}

/// Generic ellipse-vs-implicit-surface intersection via Newton refinement.
///
/// Analogous to [`super::curve_surface::circle_vs_implicit_surface`] but
/// parameterised through the ellipse's major/minor radii.
fn ellipse_vs_implicit_surface(
    ellipse: &Ellipse3,
    t_range: [f64; 2],
    f: impl Fn(DVec3) -> f64,
    geom_tol: f64,
) -> Vec<EllipseSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    use std::f64::consts::TAU;

    let cn = ellipse.normal.normalize();
    let cu = ellipse.major_dir.normalize();
    let cv = cn.cross(cu).normalize();
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;

    let pt = |theta: f64| -> DVec3 {
        ellipse.center + a * theta.cos() * cu + b * theta.sin() * cv
    };

    const N_SEEDS: usize = 64;
    let [t0, t1] = t_range;
    let span = (t1 - t0).min(TAU);

    // Sign-change detection over coarse grid
    let mut seeds: Vec<f64> = Vec::new();
    let mut prev_val = f(pt(t0));
    for i in 1..=N_SEEDS {
        let theta = t0 + span * i as f64 / N_SEEDS as f64;
        let val = f(pt(theta));
        if prev_val * val <= 0.0 {
            seeds.push(theta - span * 0.5 / N_SEEDS as f64);
        }
        prev_val = val;
    }

    // Newton refinement
    let mut hits: Vec<EllipseSurfaceHit> = Vec::new();
    const MAX_ITER: usize = 20;
    const H: f64 = TOLERANCE_ABS;
    for seed in seeds {
        let mut theta = seed;
        for _ in 0..MAX_ITER {
            let fv = f(pt(theta));
            let dfdtheta = (f(pt(theta + H)) - f(pt(theta - H))) / (2.0 * H);
            if dfdtheta.abs() < TOLERANCE_LEN_SQ_DIV_SAFE {
                break;
            }
            let delta = -fv / dfdtheta;
            theta += delta;
            if delta.abs() < eps * 0.01 {
                break;
            }
        }

        if theta < t0 - eps || theta > t1 + eps {
            continue;
        }
        let point = pt(theta);
        if f(point).abs() > eps * 10.0 {
            continue;
        }

        let duplicate = hits.iter().any(|h: &EllipseSurfaceHit| {
            (h.ellipse_param - theta).abs() < eps * 5.0
        });
        if !duplicate {
            hits.push(EllipseSurfaceHit {
                point,
                ellipse_param: theta,
            });
        }
    }
    hits
}


