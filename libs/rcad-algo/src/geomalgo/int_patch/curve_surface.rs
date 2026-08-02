// OCCT IntCurveSurface_Inter / IntAna exact analytic curve-surface intersection.
// Ported from the legacy rcad-algorithms implementation (curve_surface.rs),
// which mirrors OCCT IntAna_IntConicQuad for analytic curve/quadric pairs.
use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, CylindricalSurface, Line3, Plane, SphericalSurface};
use rcad_kernel::precision::CONFUSION as TOLERANCE_ABS;

const TOLERANCE_LEN_SQ_DIV_SAFE: f64 = 1e-30;
const TOLERANCE_TOL_SCALE_MICRO: f64 = 1e-6;

/// Hit from a curve-surface intersection.
pub struct CurveSurfaceHit {
    pub point: DVec3,
    /// Parametric value on the curve.
    pub curve_param: f64,
}

/// Intersect a line with a cylindrical surface (infinite cylinder).
/// `geom_tol` widens the accepted edge parameter interval (minimum TOLERANCE_ABS).
/// Quadratic degeneracy and discriminant tests stay strict.
pub fn intersect_line_cylinder_with_tol(
    line: &Line3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let param_eps = geom_tol.max(TOLERANCE_ABS);
    // Project to 2D perpendicular to cylinder axis
    let oc = line.origin - cyl.origin;

    // Component of direction perpendicular to axis
    let d_perp = line.direction - cyl.axis * line.direction.dot(cyl.axis);
    let oc_perp = oc - cyl.axis * oc.dot(cyl.axis);

    let a = d_perp.dot(d_perp);
    let b = 2.0 * oc_perp.dot(d_perp);
    let c = oc_perp.dot(oc_perp) - cyl.radius * cyl.radius;

    solve_quadratic_hits(a, b, c, line, t_range, param_eps)
}

/// Intersect a line with a sphere.
pub fn intersect_line_sphere_with_tol(
    line: &Line3,
    t_range: [f64; 2],
    sphere: &SphericalSurface,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let param_eps = geom_tol.max(TOLERANCE_ABS);
    let oc = line.origin - sphere.center;
    let a = line.direction.dot(line.direction);
    let b = 2.0 * oc.dot(line.direction);
    let c = oc.length_squared() - sphere.radius * sphere.radius;

    solve_quadratic_hits(a, b, c, line, t_range, param_eps)
}

fn solve_quadratic_hits(
    a: f64,
    b: f64,
    c: f64,
    line: &Line3,
    t_range: [f64; 2],
    param_eps: f64,
) -> Vec<CurveSurfaceHit> {
    let pe = param_eps.max(TOLERANCE_ABS);
    let mut hits = Vec::new();

    if a.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // Linear
        if b.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
            return hits;
        }
        let t = -c / b;
        if t >= t_range[0] - pe && t <= t_range[1] + pe {
            hits.push(CurveSurfaceHit {
                point: line.origin + line.direction * t,
                curve_param: t,
            });
        }
        return hits;
    }

    let discriminant = b * b - 4.0 * a * c;
    if discriminant < -TOLERANCE_ABS {
        return hits;
    }

    if discriminant.abs() < TOLERANCE_ABS {
        let t = -b / (2.0 * a);
        if t >= t_range[0] - pe && t <= t_range[1] + pe {
            hits.push(CurveSurfaceHit {
                point: line.origin + line.direction * t,
                curve_param: t,
            });
        }
    } else {
        let sqrt_d = discriminant.sqrt();
        for t in [(-b - sqrt_d) / (2.0 * a), (-b + sqrt_d) / (2.0 * a)] {
            if t >= t_range[0] - pe && t <= t_range[1] + pe {
                hits.push(CurveSurfaceHit {
                    point: line.origin + line.direction * t,
                    curve_param: t,
                });
            }
        }
    }
    hits
}

/// Intersect a line with a conical surface (infinite cone).
/// Nappe / parameter margins from `geom_tol`.
pub fn intersect_line_cone_with_tol(
    line: &Line3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let param_eps = geom_tol.max(TOLERANCE_ABS);
    let apex = cone.apex_point();
    let axis = cone.axis_dir();
    let co = line.origin - apex;
    let cos2 = cone.half_angle_rad.cos().powi(2);

    let d_dot_a = line.direction.dot(axis);
    let co_dot_a = co.dot(axis);

    // Point P on cone satisfies:
    //   ((P-apex)·axis)² = cos²(half_angle) * |P-apex|²
    // Substituting P = O + t*D:
    let d_d = line.direction.dot(line.direction);
    let co_d = co.dot(line.direction);
    let co_co = co.dot(co);

    let a = d_dot_a * d_dot_a - cos2 * d_d;
    let b = 2.0 * (d_dot_a * co_dot_a - cos2 * co_d);
    let c = co_dot_a * co_dot_a - cos2 * co_co;

    // Both nappes are kept: the cone's actual face (frustum) may lie on either
    // side of the apex depending on orientation, and the face's UV domain /
    // point-in-face check filters the correct one (OCCT IntAna_IntConicQuad
    // does not filter by nappe here).
    solve_quadratic_hits(a, b, c, line, t_range, param_eps)
}

/// Intersect a circle arc with a plane. Returns 0, 1, or 2 hits.
pub fn intersect_circle_plane_with_tol(
    circle: &Circle3,
    t_range: [f64; 2],
    plane: &Plane,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    // Circle parametric: P(θ) = center + radius*(u*cos(θ) + v*sin(θ)).
    // Use the circle's own local frame so the returned curve_param θ is in the
    // circle's parameter space (the range manager operates in that space).
    let u = circle.x_dir.normalize();
    let v = circle.y_dir.normalize();

    // Plane equation: (P - plane.origin) · plane.normal = 0
    let d = (circle.center - plane.origin).dot(plane.normal);
    let a_coeff = circle.radius * u.dot(plane.normal);
    let b_coeff = circle.radius * v.dot(plane.normal);

    // d + a_coeff*cos(θ) + b_coeff*sin(θ) = 0
    let r_amp = (a_coeff * a_coeff + b_coeff * b_coeff).sqrt();
    if r_amp < TOLERANCE_ABS {
        return vec![]; // circle parallel to plane
    }

    let ratio = -d / r_amp;
    if ratio.abs() > 1.0 + eps {
        return vec![];
    }
    let ratio = ratio.clamp(-1.0, 1.0);

    let phi = b_coeff.atan2(a_coeff);
    let alpha = ratio.acos();

    // The circle's parameterization may be offset (e.g. the edge range starts
    // at 3*PI/2). The raw theta is periodic with period 2π, so shift it by
    // multiples of 2π until it falls inside [t_range[0], t_range[0] + 2π),
    // then accept if it also lies within [t_range[0], t_range[1]].
    let period = 2.0 * std::f64::consts::PI;
    let mut hits = Vec::new();
    for theta_raw in [phi + alpha, phi - alpha] {
        let mut theta = theta_raw + period * ((t_range[0] - theta_raw) / period).ceil();
        if theta >= t_range[0] - eps && theta <= t_range[1] + eps {
            let point = circle.center
                + u * (circle.radius * theta.cos())
                + v * (circle.radius * theta.sin());
            hits.push(CurveSurfaceHit {
                point,
                curve_param: theta,
            });
        }
    }
    hits
}

/// Intersect a circle arc with a cylindrical surface. Returns 0-4 hits.
/// Newton refinement seeded from a coarse angle grid.
pub fn intersect_circle_cylinder_with_tol(
    circle: &Circle3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    circle_vs_implicit_surface(
        circle,
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

/// Intersect a circle arc with a spherical surface. Returns 0-2 hits.
pub fn intersect_circle_sphere_with_tol(
    circle: &Circle3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    circle_vs_implicit_surface(
        circle,
        t_range,
        |p: DVec3| -> f64 { (p - sph.center).length_squared() - sph.radius * sph.radius },
        eps,
    )
}

/// Intersect a circle arc with a conical surface. Returns 0-4 hits.
pub fn intersect_circle_cone_with_tol(
    circle: &Circle3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let cos2 = cone.half_angle_rad.cos().powi(2);
    let apex = cone.apex_point();
    let axis = cone.axis_dir();
    let eps = geom_tol.max(TOLERANCE_ABS);
    circle_vs_implicit_surface(
        circle,
        t_range,
        |p: DVec3| -> f64 {
            let v = p - apex;
            let along = v.dot(axis);
            let along2 = along * along;
            let len2 = v.length_squared();
            // Cone implicit: (v·axis)² = cos²(half) * |v|²
            along2 - cos2 * len2
        },
        eps,
    )
}

/// Generic circle-vs-implicit-surface intersection via Newton refinement.
fn circle_vs_implicit_surface(
    circle: &Circle3,
    t_range: [f64; 2],
    f: impl Fn(DVec3) -> f64,
    geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    use std::f64::consts::TAU;

    // Build a local orthonormal frame for the circle — use the circle's own
    // local frame so the returned curve_param θ is in the circle's parameter
    // space (the range manager operates in that space).
    let cu = circle.x_dir.normalize();
    let cv = circle.y_dir.normalize();

    let pt = |theta: f64| -> DVec3 {
        circle.center + circle.radius * (theta.cos() * cu + theta.sin() * cv)
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
            // Sign change — midpoint as seed
            seeds.push(theta - span * 0.5 / N_SEEDS as f64);
        }
        prev_val = val;
    }

    // Newton refinement
    let mut hits: Vec<CurveSurfaceHit> = Vec::new();
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

        // Validate within t_range and on the surface
        if theta < t0 - eps || theta > t1 + eps {
            continue;
        }
        let point = pt(theta);
        if f(point).abs() > eps * 10.0 {
            continue;
        }

        // Deduplicate
        let duplicate = hits
            .iter()
            .any(|h: &CurveSurfaceHit| (h.curve_param - theta).abs() < eps * 5.0);
        if !duplicate {
            hits.push(CurveSurfaceHit {
                point,
                curve_param: theta,
            });
        }
    }
    hits
}
