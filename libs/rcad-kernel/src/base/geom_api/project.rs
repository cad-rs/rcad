//! Closest-point projection from a 3D point onto a curve or surface.
//!
//! Analogous to OCCT `GeomAPI_ProjectPointOnCurve` and
//! `GeomAPI_ProjectPointOnSurf`.
//!
//! Low-level extrema algorithms (Extrema_ExtPElC, Extrema_ExtPC, Extrema_ExtPS)
//! live in [`crate::math::extrema`]; this module provides the GeomAPI-level
//! wrappers and analytic surface-specific projections.

use glam::{DVec2, DVec3};
use crate::geom::{Curve2d, Curve3, CurveEval, Surface3, SurfaceEval};
// Re-export the low-level extrema functions so base::geom_api::mod.rs and
// math::projection (which does `pub use ...::*`) can reach them.
pub use crate::base::extrema::{
    CurveProjection, SurfaceProjection,
    closest_point_on_curve,
    closest_point_on_surface_near,
};
// numeric_surface_projection is crate-internal (used by closest_point_on_surface)
use crate::base::extrema::numeric_surface_projection;

// ─────────────────────────────────────────────────────────────────────────────
// Curve projection (range-restricted)
// ─────────────────────────────────────────────────────────────────────────────

/// OCCT-aligned: closest_point_on_curve with range restriction [t0, t1].
///
/// Calls [`closest_point_on_curve`] then clamps parameter to [t0, t1] and
/// re-evaluates at the endpoint if the original projection fell outside range.
/// Matches OCCT GeomAPI_ProjectPointOnCurve::Init(curve, f, l).
pub fn closest_point_on_curve_range(
    curve: &Curve3,
    query: DVec3,
    t0: f64,
    t1: f64,
    n_samples: usize,
) -> CurveProjection {
    let mut result = closest_point_on_curve(curve, query, n_samples);
    if result.param < t0 || result.param > t1 {
        // OCCT GeomAPI_ProjectPointOnCurve::Init(curve, f, l) projects onto the
        // trimmed periodic curve: the parameter is wrapped by the period to stay
        // inside [f, l] (Extrema_ExtPC on a circle/ellipse considers u +/- k*Tau).
        // Clamping to an endpoint would be wrong for e.g. an arc [3*PI/2, 5*PI/2]
        // whose interior point projects to a raw parameter near 0.
        let period = match curve {
            Curve3::Circle(_) | Curve3::Ellipse(_) => Some(std::f64::consts::TAU),
            _ => None,
        };
        if let Some(p) = period {
            let k = ((t0 - result.param) / p).floor();
            for offset in [k, k + 1.0] {
                let u = result.param + offset * p;
                if u >= t0 - 1e-12 && u <= t1 + 1e-12 {
                    let pt = curve.point_at(u);
                    return CurveProjection {
                        point: pt,
                        param: u,
                        distance: (pt - query).length(),
                    };
                }
            }
        }
        let (t_min, t_max) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        let t_clamped = if result.param < t_min { t_min } else { t_max };
        let pt = curve.point_at(t_clamped);
        result = CurveProjection {
            point: pt,
            param: t_clamped,
            distance: (pt - query).length(),
        };
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface projection (analytic dispatch + fallback to numeric)
// ─────────────────────────────────────────────────────────────────────────────

/// Project the point `query` onto `surface`, returning the nearest point on the
/// surface, its (u, v) parameters, and the Euclidean distance.
///
/// Analytic surfaces (Plane, Sphere, Cylinder, Cone, Torus) use closed-form
/// formulae.  All other surfaces fall back to numerical sampling + Newton
/// refinement.
pub fn closest_point_on_surface(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    use crate::geom::*;

    match surface {
        // ── Analytic closed-form projections ──────────────────────────────────
        Surface3::Plane(plane) => {
            let d = (query - plane.origin).dot(plane.normal);
            let point = query - plane.normal * d;
            let distance = d.abs();
            let u_axis = plane.u_dir;
            let v_axis = plane.v_dir;
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
                sph.center + sph.radius * DVec3::X
            } else {
                sph.center + v / len * sph.radius
            };
            let w = (point - sph.center).normalize_or_zero();
            let u_axis = sph.ref_dir;
            let v_axis = sph.axis.cross(u_axis);
            // OCCT ElSLib::SphereD0: the axis component is R*sin(v), so
            // v = asin(w.axis) (0 = equator, +pi/2 = axis pole).
            let theta = w.dot(sph.axis).clamp(-1.0, 1.0).asin();
            let mut phi = w.dot(v_axis).atan2(w.dot(u_axis));
            // OCCT ElSLib::SphereParameters (ElSLib.cxx L1643): U is
            // normalized into [0, 2*PI) (normalizeAngle).  Without it a point
            // on the -Y side projects to u = -PI/2 and the pcurve_2d periodic
            // shift folds the arc into [2*PI, ...], so the WireSplitter's
            // closed-vertex 2D distance check compares u=2*PI (seam pcurve2)
            // against u=0 and filters the candidate (bopfuse_simple ZH6:
            // sphere equator section pcurves).
            if phi < 0.0 {
                phi += std::f64::consts::TAU;
            }
            SurfaceProjection {
                point,
                params: (phi, theta),
                distance: (point - query).length(),
            }
        }

        Surface3::Cylinder(cyl) => {
            let v = query - cyl.origin;
            let along = v.dot(cyl.axis);
            let radial = v - cyl.axis * along;
            let radial_len = radial.length();
            let point = if radial_len < 1e-14 {
                cyl.origin + cyl.axis * along + cyl.radius * cyl.ref_dir
            } else {
                cyl.origin + cyl.axis * along + radial / radial_len * cyl.radius
            };
            let u_axis = cyl.ref_dir.normalize();
            let v_axis = cyl.axis.cross(u_axis);
            // OCCT ElSLib::CylinderParameters: atan2(P·YDirection, P·XDirection)
            // on the raw offset vector (same rationale as the Cone branch).
            let theta = v.dot(v_axis).atan2(v.dot(u_axis));
            let theta = if theta < 0.0 { theta + std::f64::consts::TAU } else { theta };
            SurfaceProjection {
                point,
                params: (theta, along),
                distance: (point - query).length(),
            }
        }

        Surface3::Cone(cone) => {
            let axis = cone.axis_dir();
            // OCCT gp_Cone u=0 is the XDirection of its gp_Ax3 (the cone's
            // reference direction), preserved through rotation. any_perpendicular
            // here would break the UV mapping of rotated cones (P014).
            let x_axis = cone.ref_dir.normalize_or_zero();
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let local = query - cone.apex;
            let along = local.dot(axis);
            let radial = local - axis * along;
            let radial_len = radial.length();
            let half = cone.half_angle_rad;
            let tan_h = half.tan();
            let axial = (along + (radial_len - cone.radius) * tan_h) / (1.0 + tan_h * tan_h);
            let r_hat = if radial_len < 1e-14 { x_axis } else { radial / radial_len };
            let point = cone.apex + axis * axial + r_hat * cone.radius_at_axial(axial);
            let slant = cone.slant_from_axial(axial);
            // OCCT ElSLib::ConeParameters (ElSLib.cxx L525-556): the azimuth is
            // atan2(P·YDirection, P·XDirection) on the raw offset vector — NOT
            // on the normalized radial. Normalizing first re-introduces rounding
            // that flips the sign of a ~0 Y component (atan2(-1e-16, X) -> 2PI),
            // which makes a vertex sitting exactly on u=0 project to u=2PI and
            // breaks the WireSplitter's periodic UV comparison (P014).
            let theta = local.dot(y_axis).atan2(local.dot(x_axis));
            // OCCT gp_Cone u in [0, 2*PI]; normalize the atan2 azimuth so the
            // projected pcurves share the cone's natural UV domain (matches the
            // Cylinder branch above).
            let theta = if theta < 0.0 { theta + std::f64::consts::TAU } else { theta };
            SurfaceProjection {
                point,
                params: (theta, slant),
                distance: (point - query).length(),
            }
        }

        Surface3::Torus(torus) => {
            let v = query - torus.center;
            let along = v.dot(torus.axis);
            let radial = v - torus.axis * along;
            let radial_len = radial.length();
            let torus_ref_dir = {
                let ref_dir = if torus.axis.x.abs() > 1.0 - 1e-12 { DVec3::Z } else { DVec3::X };
                (ref_dir - torus.axis * ref_dir.dot(torus.axis)).normalize_or_zero()
            };
            let major_dir = if radial_len < 1e-14 { torus_ref_dir } else { radial / radial_len };
            let tube_center = torus.center + major_dir * torus.major_radius;
            let w = query - tube_center;
            let w_len = w.length();
            let point = if w_len < 1e-14 {
                tube_center + major_dir * torus.minor_radius
            } else {
                tube_center + w / w_len * torus.minor_radius
            };
            let u = major_dir.dot(torus_ref_dir).atan2(major_dir.dot(torus.axis.cross(torus_ref_dir)));
            let w_dir = (point - tube_center).normalize_or_zero();
            let v_param = w_dir.dot(torus.axis).atan2(w_dir.dot(major_dir));
            SurfaceProjection {
                point,
                params: (u, v_param),
                distance: (point - query).length(),
            }
        }

        // ── Planar BSpline: use analytic plane projection ──────────────────────
        Surface3::BSpline(bsp)
            if bsp.degree_u == 1
                && bsp.degree_v == 1
                && bsp.control_points.len() >= 2
                && bsp.control_points[0].len() >= 2 =>
        {
            let p00 = bsp.control_points[0][0];
            let p10 = bsp.control_points[1][0];
            let p01 = bsp.control_points[0][1];
            let du = p10 - p00;
            let dv = p01 - p00;
            let normal = du.cross(dv).normalize_or_zero();
            if normal.length_squared() > 0.5 {
                let n = normal;
                let d = (query - p00).dot(n);
                let point = query - n * d;
                let diff = point - p00;
                let u_len2 = du.length_squared();
                let v_len2 = dv.length_squared();
                let (u, v) = if u_len2 > 1e-30 && v_len2 > 1e-30 {
                    (diff.dot(du) / u_len2, diff.dot(dv) / v_len2)
                } else {
                    (0.0, 0.0)
                };
                SurfaceProjection { point, params: (u, v), distance: d.abs() }
            } else {
                numeric_surface_projection(surface, query, n_samples)
            }
        }

        // ── Numerical fallback for parametric surfaces ─────────────────────────
        _ => numeric_surface_projection(surface, query, n_samples),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// p-curve construction (GeomAPI-level)
// ─────────────────────────────────────────────────────────────────────────────

/// OCCT-aligned: IntTools_Curve::MakePCurveOnSurface.
///
/// Samples the 3D curve at `n_samples` uniformly-spaced parameter values
/// within `t_range`, projects each onto `surface` via [`closest_point_on_surface`],
/// and fits a [`Curve2d::BSpline`] through the UV points.
///
/// Returns `None` when fewer than 2 valid projections are obtained, or
/// when the curve fit fails (coincident UV points, degenerate sample set).
pub fn make_pcurve_on_surface(
    curve: &Curve3,
    t_range: [f64; 2],
    surface: &Surface3,
    n_samples: usize,
) -> Option<Curve2d> {
    use crate::math::fit::interpolate_points_2d;
    use glam::DVec2;

    let n = n_samples.max(2);
    let [t0, t1] = t_range;
    let dt = if n > 1 { (t1 - t0) / (n - 1) as f64 } else { 0.0 };

    let mut uv_pts: Vec<DVec2> = Vec::with_capacity(n);
    for i in 0..n {
        let t = t0 + dt * i as f64;
        let pt3d = curve.point_at(t);
        let proj = closest_point_on_surface(surface, pt3d, 8);
        if proj.distance > 1e-4 { continue; }
        uv_pts.push(DVec2::new(proj.params.0, proj.params.1));
    }

    if uv_pts.len() < 2 { return None; }
    uv_pts.dedup_by(|a, b| (*a - *b).length_squared() < 1e-20);

    let bspline = interpolate_points_2d(&uv_pts).ok()?;
    Some(Curve2d::BSpline(bspline))
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
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Y));
        let q = DVec3::new(3.0, 5.0, -2.0);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!((r.point - DVec3::new(3.0, 0.0, -2.0)).length() < 1e-9);
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO, axis: DVec3::Y, radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Y),
        });
        let q = DVec3::new(5.0, 0.0, 0.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        assert!((r.point - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 3.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_cylinder() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO, axis: DVec3::Y, radius: 1.0, ref_dir: DVec3::X,
        });
        let q = DVec3::new(3.0, 2.0, 0.0);
        let r = closest_point_on_surface(&cyl, q, 16);
        assert!((r.point - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 2.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_torus() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO, axis: DVec3::Y, major_radius: 3.0, minor_radius: 1.0,
        });
        let q = DVec3::new(10.0, 0.0, 0.0);
        let r = closest_point_on_surface(&torus, q, 16);
        assert!((r.point - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn project_onto_cone_returns_theta_and_slant_params() {
        let cone = Surface3::Cone(ConicalSurface::new(
            DVec3::ZERO, DVec3::Z, 2.0, 30.0_f64.to_radians(),
        ));
        let expected_slant = 4.0;
        let on_surface = match &cone {
            Surface3::Cone(surface) => surface.point_at(0.0, expected_slant),
            _ => unreachable!(),
        };
        let query_normal = match &cone {
            Surface3::Cone(surface) => surface.normal_at(0.0, expected_slant),
            _ => unreachable!(),
        };
        let q = on_surface + query_normal * 0.25;
        let r = closest_point_on_surface(&cone, q, 16);
        assert!((r.point - on_surface).length() < 5e-3);
        assert!((r.params.1 - expected_slant).abs() < 5e-3);
        let lifted = match &cone {
            Surface3::Cone(surface) => surface.point_at(r.params.0, r.params.1),
            _ => unreachable!(),
        };
        assert!((lifted - r.point).length() < 1e-6);
    }

    #[test]
    fn project_onto_bspline_surface() {
        use crate::geom::BSplineSurface;
        let surf = Surface3::BSpline(BSplineSurface {
            degree_u: 1, degree_v: 1,
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
        assert!((r.point - DVec3::new(0.5, 0.5, 0.0)).length() < 1e-4);
    }

    #[test]
    fn project_onto_cone_surface() {
        let cone = Surface3::Cone(ConicalSurface::new(
            DVec3::ZERO, DVec3::Z, 1.0, std::f64::consts::FRAC_PI_6,
        ));
        let q = DVec3::new(1.0, 0.0, 1.0);
        let r = closest_point_on_surface(&cone, q, 16);
        assert!(r.distance < 0.5);

        let q2 = DVec3::new(0.0, 0.0, 5.0);
        let r2 = closest_point_on_surface(&cone, q2, 16);
        assert!(r2.distance > 0.0);
    }

    #[test]
    fn project_onto_torus_surface() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0,
        });
        let q = DVec3::new(0.0, 0.0, 0.0);
        let r = closest_point_on_surface(&torus, q, 16);
        assert!(r.distance > 0.0);

        let q2 = DVec3::new(4.0, 0.0, 0.0);
        let r2 = closest_point_on_surface(&torus, q2, 16);
        assert!((r2.distance - 0.0).abs() < 0.1);
    }

    #[test]
    fn project_onto_cylinder_interior() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO, axis: DVec3::Z, radius: 2.0, ref_dir: DVec3::X,
        });
        let q = DVec3::new(0.0, 0.0, 1.0);
        let r = closest_point_on_surface(&cyl, q, 16);
        assert!((r.distance - 2.0).abs() < 1e-6);
    }

    #[test]
    fn project_onto_sphere_interior() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO, axis: DVec3::Z, radius: 3.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let q = DVec3::new(1.0, 1.0, 1.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        assert!(r.distance < 3.0);
    }

    #[test]
    fn project_onto_plane_offset() {
        let plane = Surface3::Plane(Plane::new(DVec3::new(0.0, 5.0, 0.0), DVec3::Y));
        let q = DVec3::new(3.0, 0.0, -2.0);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!((r.point.y - 5.0).abs() < 1e-9);
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_distant_point_onto_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let q = DVec3::new(1000.0, 1000.0, 1000.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        let expected_dist = q.length() - 1.0;
        assert!((r.distance - expected_dist).abs() < 1.0);
    }

    #[test]
    fn project_near_surface_boundary() {
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let q = DVec3::new(1.0, 2.0, 1e-10);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!(r.distance < 1e-9);
    }
}
