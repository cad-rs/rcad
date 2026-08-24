// OCCT IntCurveSurface_Inter / IntAna exact analytic curve-surface intersection.
// Ported from the legacy rcad-algorithms implementation (curve_surface.rs),
// which mirrors OCCT IntAna_IntConicQuad for analytic curve/quadric pairs.
use glam::DVec3;
use rcad_kernel::geom::{
    any_perpendicular, Circle3, ConicalSurface, CylindricalSurface, Line3, Plane,
    SphericalSurface, ToroidalSurface,
};
use rcad_kernel::math::direct_polynomial_roots::DirectPolynomialRoots;
use rcad_kernel::precision::CONFUSION as TOLERANCE_ABS;
use rcad_kernel::SurfaceEval;

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

/// Intersect a line with a torus.  OCCT IntAna_IntLinTorus
/// (IntAna_IntLinTorus.cxx L43-132): the line is re-parameterized so its
/// origin is the closest point to the torus center, transformed into the
/// torus frame, and the quartic (a4,a3,a2,a1,a0) is solved by
/// math_DirectPolynomialRoots.  Each root t is validated by re-evaluating
/// the torus at the surface parameters of the solution point
/// (ElSLib::Parameters) and checking |PSolT - PSolL|^2 <= 1e-10.
pub fn intersect_line_torus_with_tol(
    line: &Line3,
    t_range: [f64; 2],
    tor: &ToroidalSurface,
    _geom_tol: f64,
) -> Vec<CurveSurfaceHit> {
    let pl = line.origin;
    let dl = line.direction;
    let tor_loc = tor.center;
    // Reparametrize the line: set its location as nearest to the torus location.
    // OCCT: ParamOfNewPL = gp_Vec(PL, TorLoc).Dot(gp_Vec(DL)) = (TorLoc - PL)·DL.
    let param_of_new_pl = (tor_loc - pl).dot(dl);
    let new_pl = pl + param_of_new_pl * dl;
    // Transform into the torus reference frame: X/Y in the equatorial plane
    // (the torus's gp_Ax3 XDirection/YDirection — the stored ref_dir,
    // preserved through rotation), Z along the axis.
    // OCCT L57-60: NewPL.Transform(trsf) with trsf.SetTransformation(T.Position())
    // maps the world point into the torus LOCAL frame (torus center at origin).
    let z = tor.axis.normalize_or_zero();
    let x = tor.ref_dir.normalize_or_zero();
    let y = z.cross(x).normalize_or_zero();
    let rot = |v: DVec3| DVec3::new(v.dot(x), v.dot(y), v.dot(z));
    let x0 = rot(new_pl - tor_loc);
    let x1 = rot(dl);
    let (r, r2) = (tor.major_radius, tor.major_radius * tor.major_radius);
    let (r_min, r_min2) = (tor.minor_radius, tor.minor_radius * tor.minor_radius);
    let a = x1.x * x1.x + x1.y * x1.y + x1.z * x1.z;
    let b = 2.0 * (x1.x * x0.x + x1.y * x0.y + x1.z * x0.z);
    let c = x0.x * x0.x + x0.y * x0.y + x0.z * x0.z - (r2 + r_min2);
    let a4 = a * a;
    let a3 = 2.0 * a * b;
    let a2 = 2.0 * a * c + 4.0 * r2 * x1.z * x1.z + b * b;
    let a1 = 2.0 * b * c + 8.0 * r2 * x1.z * x0.z;
    let a0 = c * c + 4.0 * r2 * (x0.z * x0.z - r_min2);
    let mdpr = DirectPolynomialRoots::new_quartic(a4, a3, a2, a1, a0);
    if std::env::var("RCAD_LT_DEBUG").is_ok() {
        eprintln!("[LT] a4={:.6e} a3={:.6e} a2={:.6e} a1={:.6e} a0={:.6e} done={} n={} pl={:?} dl={:?} tor={:?} param_new={:.6} new_pl={:?}",
            a4, a3, a2, a1, a0, mdpr.is_done(), if mdpr.is_done() { mdpr.nb_solutions() } else { 0 },
            pl, dl, tor.center, param_of_new_pl, new_pl);
    }
    if !mdpr.is_done() {
        return Vec::new();
    }
    let param_eps = TOLERANCE_ABS;
    let mut hits = Vec::new();
    let n = mdpr.nb_solutions();
    let mut n_bad = 0usize;
    // OCCT mdpr.Value(i) is 1-based (math_DirectPolynomialRoots.hxx).
    for i in 1..=n {
        let t = mdpr.value(i) + param_of_new_pl;
        // OCCT ElCLib::Value(t, L) — the point on the ORIGINAL line.
        let p_sol_l = pl + dl * t;
        // ElSLib::Parameters(T, PSolL, u, v) — surface params of the solution
        // point; rcad world_to_uv matches the OCCT torus parameterization.
        let uv = tor.world_to_uv(p_sol_l);
        let p_sol_t = tor.point_at(uv.x, uv.y);
        let d2 = p_sol_t.distance_squared(p_sol_l);
        if std::env::var("RCAD_LT_DEBUG").is_ok() {
            eprintln!("[LT] root t={:.6} p_sol_l=({:.6},{:.6},{:.6}) uv=({:.6},{:.6}) p_sol_t=({:.6},{:.6},{:.6}) d2={:.3e}",
                t, p_sol_l.x, p_sol_l.y, p_sol_l.z, uv.x, uv.y, p_sol_t.x, p_sol_t.y, p_sol_t.z, d2);
        }
        if d2 > 1.0e-10 {
            n_bad += 1;
        } else if t >= t_range[0] - param_eps && t <= t_range[1] + param_eps {
            hits.push(CurveSurfaceHit {
                point: p_sol_l,
                curve_param: t,
            });
        }
    }
    if n > 0 && hits.is_empty() && n_bad == n {
        // OCCT: all solutions failed validation — not done.
        return Vec::new();
    }
    hits
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
