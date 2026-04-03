use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// A numerically sampled intersection curve.
#[derive(Debug, Clone)]
pub struct SampledCurve {
    pub points: Vec<DVec3>,
    pub is_closed: bool,
}

/// Evaluate the implicit function F(P) for a surface: F=0 on surface.
fn surface_implicit(surface: &Surface3, point: DVec3) -> f64 {
    match surface {
        Surface3::Plane(p) => (point - p.origin).dot(p.normal),
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let perp = v - c.axis * along;
            perp.length() - c.radius
        }
        Surface3::Sphere(s) => (point - s.center).length() - s.radius,
        Surface3::Cone(c) => {
            let v = point - c.apex;
            let along = v.dot(c.axis);
            let perp_len = (v - c.axis * along).length();
            perp_len - along * c.half_angle_rad.tan()
        }
        Surface3::Torus(t) => {
            let v = point - t.center;
            let along = v.dot(t.axis);
            let perp = v - t.axis * along;
            let perp_len = perp.length();
            let d = perp_len - t.major_radius;
            (d * d + along * along).sqrt() - t.minor_radius
        }
        Surface3::BSpline(_) => {
            // Approximate implicit via SurfaceEval normal_at (fallback)
            point.length() - 1.0
        }
        Surface3::LinearExtrusion(_) | Surface3::Revolution(_) => point.length() - 1.0,
    }
}

/// Compute the gradient ∇F at a point for a surface.
fn surface_gradient(surface: &Surface3, point: DVec3) -> DVec3 {
    match surface {
        Surface3::Plane(p) => p.normal,
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let perp = v - c.axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            perp / perp_len
        }
        Surface3::Sphere(s) => {
            let v = point - s.center;
            let len = v.length();
            if len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            v / len
        }
        Surface3::Cone(c) => {
            let v = point - c.apex;
            let along = v.dot(c.axis);
            let perp = v - c.axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            let tan_a = c.half_angle_rad.tan();
            perp / perp_len - c.axis * tan_a
        }
        Surface3::Torus(t) => {
            let v = point - t.center;
            let along = v.dot(t.axis);
            let perp = v - t.axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            let tube_center = t.center + perp / perp_len * t.major_radius;
            let tv = point - tube_center;
            let tv_len = tv.length();
            if tv_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            tv / tv_len
        }
        Surface3::BSpline(s) => {
            // Use SurfaceEval normal_at at (0,0) as a rough approximation
            s.normal_at(0.0, 0.0)
        }
        Surface3::LinearExtrusion(s) => s.normal_at(0.0, 0.0),
        Surface3::Revolution(s) => s.normal_at(0.0, 0.0),
    }
}

/// Project a point onto a surface using Newton iteration.
fn project_onto_surface(surface: &Surface3, point: DVec3, max_iter: usize) -> DVec3 {
    let mut p = point;
    for _ in 0..max_iter {
        let f = surface_implicit(surface, p);
        if f.abs() < TOLERANCE_ABS {
            break;
        }
        let g = surface_gradient(surface, p);
        let g_len_sq = g.length_squared();
        if g_len_sq < TOLERANCE_ABS * TOLERANCE_ABS {
            break;
        }
        p -= g * (f / g_len_sq);
    }
    p
}

/// Project a point onto the intersection of two surfaces.
fn project_onto_intersection(s1: &Surface3, s2: &Surface3, point: DVec3) -> DVec3 {
    let mut p = point;
    for _ in 0..20 {
        let f1 = surface_implicit(s1, p);
        let f2 = surface_implicit(s2, p);
        if f1.abs() < TOLERANCE_ABS && f2.abs() < TOLERANCE_ABS {
            break;
        }
        let g1 = surface_gradient(s1, p);
        let g2 = surface_gradient(s2, p);

        // Solve 2x2 system: move by λ1*g1 + λ2*g2 to satisfy both constraints
        let a11 = g1.dot(g1);
        let a12 = g1.dot(g2);
        let a22 = g2.dot(g2);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
            // Degenerate — just project onto each surface alternately
            p = project_onto_surface(s1, p, 5);
            p = project_onto_surface(s2, p, 5);
            continue;
        }
        let lambda1 = (a22 * f1 - a12 * f2) / det;
        let lambda2 = (a11 * f2 - a12 * f1) / det;
        p -= g1 * lambda1 + g2 * lambda2;
    }
    p
}

/// Find seed points for intersection curve marching by sampling one surface.
pub fn find_seed_points(
    s1: &Surface3,
    s2: &Surface3,
    sample_points: &[DVec3],
) -> Vec<DVec3> {
    let mut seeds = Vec::new();

    // Look for sign changes of F2 along the sample points
    let values: Vec<f64> = sample_points
        .iter()
        .map(|&p| surface_implicit(s2, p))
        .collect();

    for i in 0..values.len().saturating_sub(1) {
        if values[i] * values[i + 1] < 0.0 {
            // Sign change — interpolate
            let t = values[i] / (values[i] - values[i + 1]);
            let p = sample_points[i] + (sample_points[i + 1] - sample_points[i]) * t;
            let seed = project_onto_intersection(s1, s2, p);
            seeds.push(seed);
        }
    }

    seeds
}

/// March an intersection curve starting from a seed point.
/// Traces in both directions along the curve until it returns to start
/// (closed) or exits bounds.
pub fn march_intersection(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    step_size: f64,
    max_steps: usize,
    bounds_check: impl Fn(DVec3) -> bool,
) -> SampledCurve {
    let mut forward_points = march_one_direction(s1, s2, seed, step_size, max_steps, &bounds_check);
    let backward_points = march_one_direction(s1, s2, seed, -step_size, max_steps, &bounds_check);

    // Combine: reverse backward, then append forward (excluding duplicate seed)
    let mut points: Vec<DVec3> = backward_points.into_iter().rev().collect();
    if !forward_points.is_empty() {
        // Skip the seed which is duplicated
        forward_points.remove(0);
        points.extend(forward_points);
    }

    // Check closure
    let is_closed = points.len() > 2
        && points_coincide(*points.first().unwrap(), *points.last().unwrap());
    if is_closed {
        points.pop();
    }

    SampledCurve { points, is_closed }
}

fn march_one_direction(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    step_size: f64,
    max_steps: usize,
    bounds_check: &impl Fn(DVec3) -> bool,
) -> Vec<DVec3> {
    let mut points = vec![seed];
    let mut current = seed;

    for _ in 0..max_steps {
        let g1 = surface_gradient(s1, current);
        let g2 = surface_gradient(s2, current);
        let tangent = g1.cross(g2);
        let t_len = tangent.length();
        if t_len < TOLERANCE_ABS {
            break; // tangent surfaces — can't march
        }
        let dir = tangent / t_len;

        let next_raw = current + dir * step_size;
        let next = project_onto_intersection(s1, s2, next_raw);

        if !bounds_check(next) {
            break;
        }

        // Check if we've returned to start (closed curve)
        if points.len() > 3 && points_coincide(next, points[0]) {
            points.push(points[0]);
            break;
        }

        points.push(next);
        current = next;
    }

    points
}

/// Generate sample points on a cylinder surface for seed finding.
pub fn sample_cylinder(cyl: &CylindricalSurface, height_range: [f64; 2], n_theta: usize, n_h: usize) -> Vec<DVec3> {
    let u = if cyl.axis.x.abs() < 0.9 {
        cyl.axis.cross(DVec3::X).normalize()
    } else {
        cyl.axis.cross(DVec3::Y).normalize()
    };
    let v = cyl.axis.cross(u);

    let mut points = Vec::with_capacity(n_theta * n_h);
    for ih in 0..n_h {
        let h = height_range[0] + (height_range[1] - height_range[0]) * ih as f64 / (n_h - 1).max(1) as f64;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = cyl.origin + cyl.axis * h + (u * theta.cos() + v * theta.sin()) * cyl.radius;
            points.push(p);
        }
    }
    points
}

/// Generate sample points on a sphere surface for seed finding.
pub fn sample_sphere(sphere: &SphericalSurface, n_theta: usize, n_phi: usize) -> Vec<DVec3> {
    let u = if sphere.axis.x.abs() < 0.9 {
        sphere.axis.cross(DVec3::X).normalize()
    } else {
        sphere.axis.cross(DVec3::Y).normalize()
    };
    let v = sphere.axis.cross(u);

    let mut points = Vec::with_capacity(n_theta * n_phi);
    for ip in 0..n_phi {
        let phi = std::f64::consts::PI * ip as f64 / (n_phi - 1).max(1) as f64;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = sphere.center
                + sphere.radius * (sphere.axis * phi.cos() + (u * theta.cos() + v * theta.sin()) * phi.sin());
            points.push(p);
        }
    }
    points
}

/// Generate sample points on a torus surface for seed finding.
pub fn sample_torus(torus: &ToroidalSurface, n_u: usize, n_v: usize) -> Vec<DVec3> {
    let u_dir = if torus.axis.x.abs() < 0.9 {
        torus.axis.cross(DVec3::X).normalize()
    } else {
        torus.axis.cross(DVec3::Y).normalize()
    };
    let v_dir = torus.axis.cross(u_dir);

    let mut points = Vec::with_capacity(n_u * n_v);
    for iu in 0..n_u {
        let u = 2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
        let cu = u.cos();
        let su = u.sin();
        let ring_center = torus.center + (u_dir * cu + v_dir * su) * torus.major_radius;
        let ring_outward = u_dir * cu + v_dir * su;

        for iv in 0..n_v {
            let v = 2.0 * std::f64::consts::PI * iv as f64 / n_v as f64;
            let p = ring_center + (ring_outward * v.cos() + torus.axis * v.sin()) * torus.minor_radius;
            points.push(p);
        }
    }
    points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn implicit_plane() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        });
        assert!((surface_implicit(&plane, DVec3::ZERO)).abs() < TOLERANCE_ABS);
        assert!((surface_implicit(&plane, DVec3::new(0.0, 1.0, 0.0)) - 1.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn implicit_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        assert!((surface_implicit(&sphere, DVec3::new(2.0, 0.0, 0.0))).abs() < TOLERANCE_ABS);
        assert!((surface_implicit(&sphere, DVec3::new(1.0, 0.0, 0.0)) + 1.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn implicit_cylinder() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        });
        assert!((surface_implicit(&cyl, DVec3::new(3.0, 5.0, 0.0))).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn implicit_torus() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        });
        // Point on the outer equator: (6, 0, 0)
        assert!((surface_implicit(&torus, DVec3::new(6.0, 0.0, 0.0))).abs() < TOLERANCE_ABS);
        // Point on the inner equator: (4, 0, 0)
        assert!((surface_implicit(&torus, DVec3::new(4.0, 0.0, 0.0))).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn project_onto_sphere_test() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        let p = project_onto_surface(&sphere, DVec3::new(3.0, 0.0, 0.0), 20);
        assert!((p.length() - 2.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn march_sphere_cylinder_intersection() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::Y,
            radius: 1.5,
        });

        // Find a seed by sampling the cylinder and looking for sign changes on sphere
        let cyl_surf = match &cyl {
            Surface3::Cylinder(c) => c,
            _ => unreachable!(),
        };
        let samples = sample_cylinder(cyl_surf, [-1.0, 1.0], 32, 4);
        let seeds = find_seed_points(&cyl, &sphere, &samples);

        assert!(!seeds.is_empty(), "Should find at least one seed point");
        let seed = seeds[0];

        // Verify seed is approximately on both surfaces
        assert!(surface_implicit(&sphere, seed).abs() < 0.1,
            "seed not near sphere: F={}", surface_implicit(&sphere, seed));
        assert!(surface_implicit(&cyl, seed).abs() < 0.1,
            "seed not near cylinder: F={}", surface_implicit(&cyl, seed));

        let curve = march_intersection(&sphere, &cyl, seed, 0.1, 200, |_| true);
        assert!(curve.points.len() > 5,
            "Expected marched curve with several points, got {}", curve.points.len());

        // All points should be approximately on both surfaces
        for p in &curve.points {
            assert!(surface_implicit(&sphere, *p).abs() < 0.05,
                "point not on sphere: F={}", surface_implicit(&sphere, *p));
            assert!(surface_implicit(&cyl, *p).abs() < 0.05,
                "point not on cylinder: F={}", surface_implicit(&cyl, *p));
        }
    }

    #[test]
    fn sample_cylinder_test() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        let pts = sample_cylinder(&cyl, [0.0, 2.0], 8, 4);
        assert_eq!(pts.len(), 32);
        for p in &pts {
            let r = (p.x * p.x + p.z * p.z).sqrt();
            assert!((r - 1.0).abs() < TOLERANCE_ABS);
        }
    }
}
