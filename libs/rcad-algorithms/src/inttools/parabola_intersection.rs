//! Analytic intersection of a Parabola3 with analytic surfaces.
//!
//! OCCT reference: IntAna_IntConicQuad — BOPAlgo_PaveFiller_5.cxx parabola branch.
//! OCCT dispatches Parabola × {Plane, Cylinder, Cone, Sphere} through
//! PerformConicSurf → IntAna_IntConicQuad.
//!
//! ✅ OCCT-aligned: Parabola × Plane uses the quadratic solve that matches
//!    IntAna_IntConicQuad's implicit-form approach.
//! ⏳ Parabola × {Cylinder, Cone, Sphere} fall back to numeric marching
//!    (conic-implicit pairing is rare in boolean operations).

use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// Hit from a parabola-surface intersection.
pub struct ParabolaSurfaceHit {
    pub point: DVec3,
    /// Parametric value on the parabola.
    pub parabola_param: f64,
}

/// Intersect a parabola arc with a plane.
///
/// OCCT reference: IntAna_IntConicQuad — the parabola is a conic (Geom_Conic)
/// evaluated at parameter t. Substituting:
///
///   P(t) = vertex + (t²/(2p))·axis_dir + t·perp_dir
///
/// into the plane implicit n·(X - origin) = 0 gives:
///
///   n·(vertex - origin) + (t²/(2p))·(n·axis_dir) + t·(n·perp_dir) = 0
///
/// i.e. A·t² + B·t + C = 0  where:
///   C = n·(vertex - origin)
///   B = n·perp_dir
///   A = n·axis_dir / (2p)
///
/// Returns 0, 1, or 2 hits within `t_range`.
///
/// ✅ OCCT-aligned: same analytic method as IntAna_IntConicQuad for planes.
pub fn intersect_parabola_plane(
    parabola: &Parabola3,
    t_range: [f64; 2],
    plane: &Plane,
) -> Vec<ParabolaSurfaceHit> {
    intersect_parabola_plane_with_tol(parabola, t_range, plane, TOLERANCE_ABS)
}

/// Same as [`intersect_parabola_plane`] with parameter margin from `geom_tol`.
pub fn intersect_parabola_plane_with_tol(
    parabola: &Parabola3,
    t_range: [f64; 2],
    plane: &Plane,
    geom_tol: f64,
) -> Vec<ParabolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);

    // Local frame: axis_dir and perp_dir (normal × axis_dir)
    let axis_dir = parabola.axis_dir.normalize();
    let cn = parabola.normal.normalize();
    let perp_dir = cn.cross(axis_dir).normalize();

    // Plane implicit: (P - plane.origin)·n = 0
    // P(t) = vertex + (t²/(2p))·axis_dir + t·perp_dir
    let p = parabola.focal_param;
    if p.abs() < TOLERANCE_ABS {
        return vec![]; // degenerate parabola
    }

    let n = plane.normal;
    let c_term = (parabola.vertex - plane.origin).dot(n);
    let b_coeff = perp_dir.dot(n);
    let a_coeff = axis_dir.dot(n) / (2.0 * p);

    // A·t² + B·t + C = 0
    if a_coeff.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // Degenerate case: axis_dir is perpendicular to plane normal (A ≈ 0)
        if b_coeff.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
            return vec![]; // parabola parallel to plane
        }
        let t = -c_term / b_coeff;
        if t >= t_range[0] - eps && t <= t_range[1] + eps {
            let point = parabola.vertex + (t * t / (2.0 * p)) * axis_dir + t * perp_dir;
            return vec![ParabolaSurfaceHit {
                point,
                parabola_param: t,
            }];
        }
        return vec![];
    }

    let disc = b_coeff * b_coeff - 4.0 * a_coeff * c_term;
    if disc < -TOLERANCE_ABS {
        return vec![];
    }

    let mut hits = Vec::new();
    if disc.abs() < TOLERANCE_ABS {
        let t = -b_coeff / (2.0 * a_coeff);
        if t >= t_range[0] - eps && t <= t_range[1] + eps {
            let point = parabola.vertex + (t * t / (2.0 * p)) * axis_dir + t * perp_dir;
            hits.push(ParabolaSurfaceHit {
                point,
                parabola_param: t,
            });
        }
    } else {
        let sqrt_d = disc.sqrt();
        for t in [(-b_coeff - sqrt_d) / (2.0 * a_coeff), (-b_coeff + sqrt_d) / (2.0 * a_coeff)] {
            if t >= t_range[0] - eps && t <= t_range[1] + eps {
                let point = parabola.vertex + (t * t / (2.0 * p)) * axis_dir + t * perp_dir;
                hits.push(ParabolaSurfaceHit {
                    point,
                    parabola_param: t,
                });
            }
        }
    }
    hits
}

/// Intersect a parabola arc with a cylindrical surface.  Numeric fallback.
///
/// ⏳ Partially aligned: OCCT IntAna_IntConicQuad handles cylinder, but
/// this is rare in boolean operations — the numeric marching fallback
/// in `intersect_edge_face_numeric` handles it adequately.
pub fn intersect_parabola_cylinder(
    parabola: &Parabola3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
) -> Vec<ParabolaSurfaceHit> {
    intersect_parabola_cylinder_with_tol(parabola, t_range, cyl, TOLERANCE_ABS)
}

/// Same as [`intersect_parabola_cylinder`] with margins from `geom_tol`.
pub fn intersect_parabola_cylinder_with_tol(
    parabola: &Parabola3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
    geom_tol: f64,
) -> Vec<ParabolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    parabola_vs_implicit_surface(
        parabola,
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

/// Intersect a parabola arc with a spherical surface.  Numeric fallback.
pub fn intersect_parabola_sphere(
    parabola: &Parabola3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
) -> Vec<ParabolaSurfaceHit> {
    intersect_parabola_sphere_with_tol(parabola, t_range, sph, TOLERANCE_ABS)
}

/// Same as [`intersect_parabola_sphere`] with margins from `geom_tol`.
pub fn intersect_parabola_sphere_with_tol(
    parabola: &Parabola3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
    geom_tol: f64,
) -> Vec<ParabolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    parabola_vs_implicit_surface(
        parabola,
        t_range,
        |p: DVec3| -> f64 { (p - sph.center).length_squared() - sph.radius * sph.radius },
        eps,
    )
}

/// Intersect a parabola arc with a conical surface.  Numeric fallback.
pub fn intersect_parabola_cone(
    parabola: &Parabola3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
) -> Vec<ParabolaSurfaceHit> {
    intersect_parabola_cone_with_tol(parabola, t_range, cone, TOLERANCE_ABS)
}

/// Same as [`intersect_parabola_cone`] with margins from `geom_tol`.
pub fn intersect_parabola_cone_with_tol(
    parabola: &Parabola3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
    geom_tol: f64,
) -> Vec<ParabolaSurfaceHit> {
    let cos2 = cone.half_angle_rad.cos().powi(2);
    let apex = cone.apex_point();
    let axis = cone.axis_dir();
    let eps = geom_tol.max(TOLERANCE_ABS);
    parabola_vs_implicit_surface(
        parabola,
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

/// Generic parabola-vs-implicit-surface intersection via Newton refinement.
fn parabola_vs_implicit_surface(
    parabola: &Parabola3,
    t_range: [f64; 2],
    f: impl Fn(DVec3) -> f64,
    geom_tol: f64,
) -> Vec<ParabolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);

    let axis_dir = parabola.axis_dir.normalize();
    let cn = parabola.normal.normalize();
    let perp_dir = cn.cross(axis_dir).normalize();
    let p = parabola.focal_param;

    if p.abs() < TOLERANCE_ABS {
        return vec![];
    }

    let pt = |t: f64| -> DVec3 {
        parabola.vertex + (t * t / (2.0 * p)) * axis_dir + t * perp_dir
    };

    const N_SEEDS: usize = 64;
    let [t0, t1] = t_range;
    let span = t1 - t0;

    if span <= eps {
        return vec![];
    }

    // Sign-change detection over coarse grid
    let mut seeds: Vec<f64> = Vec::new();
    let mut prev_val = f(pt(t0));
    for i in 1..=N_SEEDS {
        let t = t0 + span * i as f64 / N_SEEDS as f64;
        let val = f(pt(t));
        if prev_val * val <= 0.0 {
            seeds.push(t - span * 0.5 / N_SEEDS as f64);
        }
        prev_val = val;
    }

    // Newton refinement
    let mut hits: Vec<ParabolaSurfaceHit> = Vec::new();
    const MAX_ITER: usize = 20;
    const H: f64 = TOLERANCE_ABS;
    for seed in seeds {
        let mut t = seed;
        for _ in 0..MAX_ITER {
            let fv = f(pt(t));
            let dfdt = (f(pt(t + H)) - f(pt(t - H))) / (2.0 * H);
            if dfdt.abs() < TOLERANCE_LEN_SQ_DIV_SAFE {
                break;
            }
            let delta = -fv / dfdt;
            t += delta;
            if delta.abs() < eps * 0.01 {
                break;
            }
        }

        if t < t0 - eps || t > t1 + eps {
            continue;
        }
        let point = pt(t);
        if f(point).abs() > eps * 10.0 {
            continue;
        }

        let duplicate = hits.iter().any(|h: &ParabolaSurfaceHit| {
            (h.parabola_param - t).abs() < eps * 5.0
        });
        if !duplicate {
            hits.push(ParabolaSurfaceHit {
                point,
                parabola_param: t,
            });
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    #[test]
    fn parabola_through_plane() {
        // Parabola: vertex at origin, axis along +X, opening right
        // P(t) = (t²/(2p), t, 0) where p = 2
        let parabola = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        };
        // Plane at X = 2: t²/4 = 2 → t² = 8 → t = ±2√2
        let plane = Plane {
            origin: DVec3::new(2.0, 0.0, 0.0),
            normal: DVec3::X,
        };
        let hits = intersect_parabola_plane(&parabola, [-5.0, 5.0], &plane);
        assert_eq!(hits.len(), 2, "parabola should intersect plane at 2 points");
        for h in &hits {
            assert!((h.point.x - 2.0).abs() < TOLERANCE_ABS);
        }
    }

    #[test]
    fn parabola_tangent_to_plane() {
        let parabola = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        };
        // Plane at X = 0 passes through vertex: t = 0 is the only solution
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };
        let hits = intersect_parabola_plane(&parabola, [-5.0, 5.0], &plane);
        // P(t) = (t²/4, t, 0). X=0 → t=0 is the only solution.
        // But the plane goes through the vertex, so B*0 + C = 0 with C=0.
        // A*t² + B*t + C = 0 → A*t² + B*t = 0 → t*(A*t + B) = 0, so t=0 and t=-B/A.
        // A = n·axis_dir/(2p) = 1/(4), B = n·perp_dir = 0, C = 0.
        // → (1/4)*t² = 0 → t=0 only (double root).
        assert_eq!(hits.len(), 1, "parabola vertex is tangent to plane");
        assert!((hits[0].parabola_param).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn parabola_misses_plane() {
        let parabola = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        };
        // Plane at X = -1: parabola only exists at X >= 0
        let plane = Plane {
            origin: DVec3::new(-1.0, 0.0, 0.0),
            normal: DVec3::X,
        };
        let hits = intersect_parabola_plane(&parabola, [-5.0, 5.0], &plane);
        assert!(hits.is_empty(), "parabola should miss plane behind vertex");
    }

    #[test]
    fn parabola_plane_linear_case() {
        // Parabola with axis perpendicular to plane normal → A ≈ 0 → linear
        let parabola = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        };
        // Plane normal = (0, 1, 0): plane is Y = 2
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };
        // Parabola: P(t) = (t²/4, t, 0). Plane: Y = 2.
        // A = n·axis_dir/(2p) = 0/(4) = 0
        // B = n·perp_dir = 1 (perp_dir = Y)
        // C = n·(vertex-origin) = -2
        // → 0*t² + 1*t - 2 = 0 → t = 2
        let hits = intersect_parabola_plane(&parabola, [-5.0, 5.0], &plane);
        assert_eq!(hits.len(), 1, "should have one linear intersection");
        assert!((hits[0].parabola_param - 2.0).abs() < TOLERANCE_ABS);
    }

    /// Parabola × cylinder — at least runs without crashing
    #[test]
    fn parabola_cylinder_no_crash() {
        let parabola = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::new(0.0, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 10.0,
        };
        let hits = intersect_parabola_cylinder(&parabola, [-5.0, 5.0], &cyl);
        // Just verify it runs
        assert!(hits.len() <= 4);
    }

    /// Parabola × sphere — at least runs without crashing
    #[test]
    fn parabola_sphere_no_crash() {
        let parabola = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        };
        let sph = SphericalSurface {
            center: DVec3::new(0.0, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 5.0,
        };
        let hits = intersect_parabola_sphere(&parabola, [-5.0, 5.0], &sph);
        assert!(hits.len() <= 2);
    }
}
