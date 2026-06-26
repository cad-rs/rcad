//! OCCT reference: IntAna_IntConicQuad only handles conic × Plane.
//! OCCT dispatches Hyperbola × {Cylinder, Cone, Sphere} through generic
//! numeric EF intersection (IntTools_EdgeFace), same as rcad.
//!
//! ✅ Hyperbola × Plane: analytic solve via cosh/sinh → quadratic in e^t.
//! ✅ Hyperbola × {Cylinder, Cone, Sphere}: numeric marching fallback,
//!    matching OCCT's generic EF path.

use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// Hit from a hyperbola-surface intersection.
pub struct HyperbolaSurfaceHit {
    pub point: DVec3,
    /// Parametric value on the hyperbola.
    pub hyperbola_param: f64,
}

/// Intersect a hyperbola arc with a plane.
///
/// OCCT reference: IntAna_IntConicQuad hyperbola branch.
/// The hyperbola is P(t) = center + a·cosh(t)·u + b·sinh(t)·v.
///
/// Substituting into n·(X - origin) = 0:
///
///   C + A·cosh(t) + B·sinh(t) = 0
///
/// where:
///   C = n·(center - origin)
///   A = a·(n·u)   (semi-major axis contribution)
///   B = b·(n·v)   (semi-minor axis contribution)
///
/// In exponential form (cosh(t) = (e^t+e^{-t})/2, sinh(t) = (e^t-e^{-t})/2):
///
///   2C + A·(e^t+e^{-t}) + B·(e^t-e^{-t}) = 0
///   2C + (A+B)·e^t + (A-B)·e^{-t} = 0
///
/// Multiply by e^t:
///   (A+B)·e^{2t} + 2C·e^t + (A-B) = 0
///
/// Let u = e^t > 0:
///   (A+B)·u² + 2C·u + (A-B) = 0
///
/// Returns 0, 1, or 2 hits within `t_range`.
///
/// ✅ OCCT-aligned: same analytic method as IntAna_IntConicQuad for planes.
pub fn intersect_hyperbola_plane(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    plane: &Plane,
) -> Vec<HyperbolaSurfaceHit> {
    intersect_hyperbola_plane_with_tol(hyperbola, t_range, plane, TOLERANCE_ABS)
}

/// Same as [`intersect_hyperbola_plane`] with parameter margin from `geom_tol`.
pub fn intersect_hyperbola_plane_with_tol(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    plane: &Plane,
    geom_tol: f64,
) -> Vec<HyperbolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);

    // Local frame
    let u = hyperbola.major_dir.normalize();
    let a = hyperbola.semi_major;
    let b = hyperbola.semi_minor;
    let cn = hyperbola.normal.normalize();
    let v = cn.cross(u).normalize();
    let n = plane.normal;

    // C + A·cosh(t) + B·sinh(t) = 0
    let c_term = (hyperbola.center - plane.origin).dot(n);
    let coeff_a = a * u.dot(n);
    let coeff_b = b * v.dot(n);

    if coeff_a.abs() < TOLERANCE_ABS && coeff_b.abs() < TOLERANCE_ABS {
        // Hyperbola plane is parallel to the intersection plane
        if c_term.abs() < TOLERANCE_ABS {
            // Hyperbola lies in the plane — infinite intersection, return nothing
            // (the caller should handle this as a special case)
            return vec![];
        }
        return vec![];
    }

    // Solve (A+B)·u² + 2C·u + (A-B) = 0, where u = e^t > 0
    let p_quad = coeff_a + coeff_b; // A+B
    let q_quad = 2.0 * c_term; // 2C
    let r_quad = coeff_a - coeff_b; // A-B

    let mut hits = Vec::new();

    // Helper to add hit if t is within range
    let mut try_add_hit = |t: f64| {
        if t.is_finite() && t >= t_range[0] - eps && t <= t_range[1] + eps {
            let point = hyperbola.center
                + u * (a * t.cosh())
                + v * (b * t.sinh());
            hits.push(HyperbolaSurfaceHit {
                point,
                hyperbola_param: t,
            });
        }
    };

    if p_quad.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // (A+B) ≈ 0: linear in u = e^t: 2C·u + (A-B) = 0
        if q_quad.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
            return vec![]; // fully degenerate
        }
        let u_val = -r_quad / q_quad;
        if u_val > eps {
            try_add_hit(u_val.ln());
        }
        return hits;
    }

    let disc = q_quad * q_quad - 4.0 * p_quad * r_quad;

    if disc < -TOLERANCE_ABS {
        return vec![];
    }

    if disc.abs() < TOLERANCE_ABS {
        let u_val = -q_quad / (2.0 * p_quad);
        if u_val > eps {
            try_add_hit(u_val.ln());
        }
    } else {
        let sqrt_d = disc.sqrt();
        for u_val in [(-q_quad - sqrt_d) / (2.0 * p_quad), (-q_quad + sqrt_d) / (2.0 * p_quad)] {
            if u_val > eps {
                try_add_hit(u_val.ln());
            }
        }
    }

    hits
}

/// Intersect a hyperbola arc with a cylindrical surface.  Numeric fallback.
///
/// ✅ Partially aligned: OCCT IntAna_IntConicQuad handles cylinder, but
/// this is extremely rare in boolean operations — the numeric marching
/// fallback in `intersect_edge_face_numeric` handles it.
pub fn intersect_hyperbola_cylinder(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
) -> Vec<HyperbolaSurfaceHit> {
    intersect_hyperbola_cylinder_with_tol(hyperbola, t_range, cyl, TOLERANCE_ABS)
}

/// Same as [`intersect_hyperbola_cylinder`] with margins from `geom_tol`.
pub fn intersect_hyperbola_cylinder_with_tol(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
    geom_tol: f64,
) -> Vec<HyperbolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    hyperbola_vs_implicit_surface(
        hyperbola,
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

/// Intersect a hyperbola arc with a spherical surface.  Numeric fallback.
pub fn intersect_hyperbola_sphere(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
) -> Vec<HyperbolaSurfaceHit> {
    intersect_hyperbola_sphere_with_tol(hyperbola, t_range, sph, TOLERANCE_ABS)
}

/// Same as [`intersect_hyperbola_sphere`] with margins from `geom_tol`.
pub fn intersect_hyperbola_sphere_with_tol(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
    geom_tol: f64,
) -> Vec<HyperbolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    hyperbola_vs_implicit_surface(
        hyperbola,
        t_range,
        |p: DVec3| -> f64 { (p - sph.center).length_squared() - sph.radius * sph.radius },
        eps,
    )
}

/// Intersect a hyperbola arc with a conical surface.  Numeric fallback.
pub fn intersect_hyperbola_cone(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
) -> Vec<HyperbolaSurfaceHit> {
    intersect_hyperbola_cone_with_tol(hyperbola, t_range, cone, TOLERANCE_ABS)
}

/// Same as [`intersect_hyperbola_cone`] with margins from `geom_tol`.
pub fn intersect_hyperbola_cone_with_tol(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
    geom_tol: f64,
) -> Vec<HyperbolaSurfaceHit> {
    let cos2 = cone.half_angle_rad.cos().powi(2);
    let apex = cone.apex_point();
    let axis = cone.axis_dir();
    let eps = geom_tol.max(TOLERANCE_ABS);
    hyperbola_vs_implicit_surface(
        hyperbola,
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

/// Generic hyperbola-vs-implicit-surface intersection via Newton refinement.
fn hyperbola_vs_implicit_surface(
    hyperbola: &Hyperbola3,
    t_range: [f64; 2],
    f: impl Fn(DVec3) -> f64,
    geom_tol: f64,
) -> Vec<HyperbolaSurfaceHit> {
    let eps = geom_tol.max(TOLERANCE_ABS);

    let cu = hyperbola.major_dir.normalize();
    let a = hyperbola.semi_major;
    let b = hyperbola.semi_minor;
    let cn = hyperbola.normal.normalize();
    let cv = cn.cross(cu).normalize();

    let pt = |t: f64| -> DVec3 {
        hyperbola.center + a * t.cosh() * cu + b * t.sinh() * cv
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
    let mut hits: Vec<HyperbolaSurfaceHit> = Vec::new();
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

        let duplicate = hits.iter().any(|h: &HyperbolaSurfaceHit| {
            (h.hyperbola_param - t).abs() < eps * 5.0
        });
        if !duplicate {
            hits.push(HyperbolaSurfaceHit {
                point,
                hyperbola_param: t,
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
    fn hyperbola_through_plane() {
        // Hyperbola: center at origin, major_dir = X, semi_major=2, semi_minor=1
        // P(t) = (2*cosh(t), 1*sinh(t), 0)
        let hyperbola = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 2.0,
            semi_minor: 1.0,
        };
        // Plane Y = 1: intersects at t where sinh(t) = 1 → t = asinh(1)
        let plane = Plane {
            origin: DVec3::new(0.0, 1.0, 0.0),
            normal: DVec3::Y,
        };
        let hits = intersect_hyperbola_plane(&hyperbola, [-5.0, 5.0], &plane);
        assert!(!hits.is_empty(), "hyperbola should intersect plane");
        for h in &hits {
            assert!((h.point.y - 1.0).abs() < TOLERANCE_ABS);
        }
    }

    #[test]
    fn hyperbola_misses_plane() {
        let hyperbola = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 2.0,
            semi_minor: 1.0,
        };
        // Plane at Y = 100 — far from the hyperbola branch
        let plane = Plane {
            origin: DVec3::new(0.0, 100.0, 0.0),
            normal: DVec3::Y,
        };
        let hits = intersect_hyperbola_plane(&hyperbola, [-5.0, 5.0], &plane);
        assert!(hits.is_empty(), "hyperbola should miss plane");
    }

    #[test]
    fn hyperbola_hits_both_branches() {
        // Hyperbola with A >> B: P(t) = (3*cosh(t), 2*sinh(t), 0)
        let hyperbola = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };
        // Plane X = 5
        // 3*cosh(t) = 5 → cosh(t) = 5/3 → t = ±acosh(5/3)
        let plane = Plane {
            origin: DVec3::new(5.0, 0.0, 0.0),
            normal: DVec3::X,
        };
        let hits = intersect_hyperbola_plane(&hyperbola, [-5.0, 5.0], &plane);
        assert_eq!(hits.len(), 2, "both hyperbola branches should intersect");
        for h in &hits {
            assert!((h.point.x - 5.0).abs() < TOLERANCE_ABS);
        }
    }

    #[test]
    fn hyperbola_tangent_to_plane() {
        // Hyperbola: P(t) = (cosh(t), 0.5*sinh(t), 0)
        let hyperbola = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 1.0,
            semi_minor: 0.5,
        };
        // Plane at X = 1: cosh(t) = 1 → t = 0 (tangent at vertex)
        // cosh(0)=1, sinh(0)=0 → P(0) = (1, 0, 0)
        let plane = Plane {
            origin: DVec3::new(1.0, 0.0, 0.0),
            normal: DVec3::X,
        };
        let hits = intersect_hyperbola_plane(&hyperbola, [-5.0, 5.0], &plane);
        assert_eq!(hits.len(), 1, "hyperbola should be tangent at vertex");
        assert!((hits[0].hyperbola_param).abs() < TOLERANCE_ABS);
        assert!((hits[0].point.x - 1.0).abs() < TOLERANCE_ABS);
    }

    /// Hyperbola × cylinder — at least runs without crashing
    #[test]
    fn hyperbola_cylinder_no_crash() {
        let hyperbola = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 1.0,
            semi_minor: 0.5,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 5.0,
        };
        let hits = intersect_hyperbola_cylinder(&hyperbola, [-5.0, 5.0], &cyl);
        assert!(hits.len() <= 4);
    }

    /// Hyperbola × sphere — at least runs without crashing
    #[test]
    fn hyperbola_sphere_no_crash() {
        let hyperbola = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 1.0,
            semi_minor: 0.5,
        };
        let sph = SphericalSurface {
            center: DVec3::new(0.0, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 3.0,
        };
        let hits = intersect_hyperbola_sphere(&hyperbola, [-5.0, 5.0], &sph);
        assert!(hits.len() <= 2);
    }
}
