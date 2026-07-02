//! Analytic derivation of 2D parametric curves (PCurves) for surface-surface
//! intersection results.
//!
//! Each function takes a 3D intersection curve together with surface geometry
//! and returns the exact [`Curve2d`] that represents that intersection in the
//! surface's (u, v) parameter domain.

use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::fit::interpolate_points_2d;
use rcad_kernel::geom::{
    Circle2d, Circle3, ConicalSurface, Curve2d, Curve2dEval, CurveEval, CylindricalSurface,
    Ellipse2d, Ellipse3, Hyperbola3, Line2d, Line3, Parabola3, Plane, SphericalSurface, Surface3, SurfaceEval,
    any_perpendicular,
};
use rcad_kernel::projection::closest_point_on_surface;

// ─────────────────────────────────────────────────────────────────────────────
// Plane functions
// ─────────────────────────────────────────────────────────────────────────────

/// Project a [`Circle3`] onto a [`Plane`]'s (u, v) domain.
///
/// Uses `any_perpendicular(plane.normal)` as the u-axis and
/// `plane.normal × u_axis` as the v-axis, matching [`Plane::point_at`].
///
/// If the circle lies in the plane (its normal is parallel to the plane
/// normal), the result is an analytic [`Circle2d`].  Otherwise the circle
/// projects to a general conic and is approximated with a [`BSplineCurve2`]
/// built from 33 sampled points.
pub fn circle_pcurve_on_plane(circle: &Circle3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis);

    // Test whether the circle lies in the plane.
    let normal_dot = circle
        .normal
        .normalize()
        .dot(plane.normal.normalize())
        .abs();
    if (normal_dot - 1.0).abs() < TOLERANCE_MESH_LEGACY {
        // Circle lies in the plane → analytic Circle2d.
        let diff = circle.center - plane.origin;
        let center_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));
        return Curve2d::Circle(Circle2d { center: center_2d, x_dir: DVec2::X, y_dir: DVec2::Y, radius: circle.radius,
         });
    }

    // Oblique case: sample the circle and project each point into the plane.
    let n_samples = 33_usize;
    let pts: Vec<DVec2> = (0..n_samples)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / (n_samples - 1) as f64;
            let p3 = circle.point_at(t);
            // Project onto the plane (drop the normal component).
            let diff = p3 - plane.origin;
            DVec2::new(diff.dot(u_axis), diff.dot(v_axis))
        })
        .collect();

    let mut bspline = interpolate_points_2d(&pts).expect("circle samples should not be degenerate");
    // OCCT-aligned: rescale knot vector from [0, 1] to [0, TAU] to match
    // the 3D circle curve's parameter range.
    let tau = std::f64::consts::TAU;
    for k in &mut bspline.knots { *k *= tau; }
    Curve2d::BSpline(bspline)
}

/// Project an [`Ellipse3`] onto a [`Plane`]'s (u, v) domain.
///
/// Returns an analytic [`Ellipse2d`] with the projected center, major
/// direction, and radii (unchanged — projection along a parallel normal
/// preserves semi-axes when the ellipse is coplanar with the plane).
pub fn ellipse_pcurve_on_plane(ellipse: &Ellipse3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis);

    let diff = ellipse.center - plane.origin;
    let center_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));

    let major_proj = DVec2::new(ellipse.major_dir.dot(u_axis), ellipse.major_dir.dot(v_axis));
    let major_dir_2d = if major_proj.length() > TOLERANCE_LEN_MIN {
        major_proj.normalize()
    } else {
        DVec2::X
    };

    Curve2d::Ellipse(Ellipse2d {
        center: center_2d,
        major_dir: major_dir_2d,
        major_radius: ellipse.major_radius,
        minor_radius: ellipse.minor_radius,
    })
}

/// Project a [`Line3`] onto a [`Plane`]'s (u, v) domain.
///
/// Returns a [`Line2d`] whose origin and direction are the projections of the
/// 3D line's origin and direction into the plane's parameter space.
pub fn line_pcurve_on_plane(line: &Line3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis);

    let diff = line.origin - plane.origin;
    let origin_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));

    let dir_2d = DVec2::new(line.direction.dot(u_axis), line.direction.dot(v_axis));
    let direction_2d = if dir_2d.length() > TOLERANCE_LEN_MIN {
        dir_2d.normalize()
    } else {
        DVec2::X
    };

    Curve2d::Line(Line2d {
        origin: origin_2d,
        direction: direction_2d,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Sphere functions
// ─────────────────────────────────────────────────────────────────────────────

/// OCCT-aligned: ProjLib_Sphere::Project(gp_Lin) — ProjLib_Sphere.cxx L181-184.
/// OCCT implementation is a stub: sets myType = GeomAbs_OtherCurve and relies on
/// the generic ProjLib_Projector::BuildResult (sampling → BSpline fit) for the
/// actual pcurve. rcad matches this by calling fallback_pcurve_by_projection.
pub fn line_pcurve_on_sphere(line: &Line3, sphere: &SphericalSurface) -> Curve2d {
    use rcad_kernel::geom::Curve3;
    let t_range = line.default_domain();
    // If domain is unbounded, clamp to a reasonable range
    let range = if t_range[0].is_finite() && t_range[1].is_finite() {
        t_range
    } else {
        [-1e3, 1e3]
    };
    fallback_pcurve_by_projection(&Curve3::Line(*line), &range, &Surface3::Sphere(*sphere))
}

/// OCCT-aligned: ProjLib_Sphere::Project(gp_Elips) — ProjLib_Sphere.cxx L186-189.
/// OCCT: stub sets myType = GeomAbs_OtherCurve, generic BuildResult handles it.
pub fn ellipse_pcurve_on_sphere(ellipse: &Ellipse3, sphere: &SphericalSurface) -> Curve2d {
    use rcad_kernel::geom::Curve3;
    let t_range = ellipse.default_domain();
    fallback_pcurve_by_projection(&Curve3::Ellipse(*ellipse), &t_range, &Surface3::Sphere(*sphere))
}

/// OCCT-aligned: ProjLib_Sphere::Project(gp_Parab) — stub → generic.
pub fn parabola_pcurve_on_sphere(parabola: &Parabola3, sphere: &SphericalSurface) -> Curve2d {
    use rcad_kernel::geom::Curve3;
    let t_range = parabola.default_domain();
    fallback_pcurve_by_projection(&Curve3::Parabola(*parabola), &t_range, &Surface3::Sphere(*sphere))
}

/// OCCT-aligned: ProjLib_Sphere::Project(gp_Hypr) — stub → generic.
pub fn hyperbola_pcurve_on_sphere(hyperbola: &Hyperbola3, sphere: &SphericalSurface) -> Curve2d {
    use rcad_kernel::geom::Curve3;
    let t_range = hyperbola.default_domain();
    fallback_pcurve_by_projection(&Curve3::Hyperbola(*hyperbola), &t_range, &Surface3::Sphere(*sphere))
}

/// ✅ OCCT-aligned: ProjLib_Sphere::Project(gp_Circ) — form-aligned pcurve.
///
/// OCCT ProJLib_Sphere_1.cxx L97-179 handles isoparametric circles analytically
/// (isIsoU/isIsoV → Line2d), wrapping the result in Geom2d_TrimmedCurve at the caller.
/// For non-isoparametric circles it falls through to general approximation.
///
/// rcad uses a unified 33-point BSpline fit with knot rescaling to [0, TAU] which
/// correctly handles all cases. The analytic Line2d shortcuts are omitted because
/// rcad's Line2d evaluation lacks the TrimmedCurve2 domain clipping that OCCT's
/// Geom2d_TrimmedCurve provides (point_at(t) = origin + t * direction over the full
/// [0, TAU] 3D parameter range would go outside sphere UV bounds for meridian lines).
///
/// BSpline fitting with domain-correct knots is equivalent for practical purposes and
/// does not panic where bare-Line2d would.
pub fn circle_pcurve_on_sphere(circle: &Circle3, sphere: &SphericalSurface) -> Curve2d {
    let u_ax = any_perpendicular(circle.normal).normalize();
    let v_ax = circle.normal.cross(u_ax).normalize();
    let n_samp = 33_usize;
    let mut pts: Vec<DVec2> = (0..n_samp)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / (n_samp - 1) as f64;
            let p3 = circle.center + circle.radius * (t.cos() * u_ax + t.sin() * v_ax);
            sphere.world_to_uv(p3)
        })
        .collect();
    // Unwrap seam discontinuities at U wrap-around
    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut pts[i..] { p.x -= std::f64::consts::TAU; }
        } else if du < -std::f64::consts::PI {
            for p in &mut pts[i..] { p.x += std::f64::consts::TAU; }
        }
    }
    match interpolate_points_2d(&pts) {
        Ok(mut bspline) => {
            if std::env::var("RCAD_DEBUG_PCURVE").is_ok() {
                eprintln!("[DBG_PCURVE] sphere pcurve: {} pts, BSpline, t_range=[{:.4},{:.4}]",
                    n_samp, bspline.knots.first().unwrap_or(&0.0), bspline.knots.last().unwrap_or(&0.0));
            }
            // OCCT-aligned: interpolate_points_2d produces chords in [0, 1],
            // but the 3D circle curve is parameterized on [0, TAU].  Rescale
            // the knot vector to [0, TAU] so point_at(t) with t from the 3D
            // curve's parameter range maps to the correct UV position.
            let tau = std::f64::consts::TAU;
            for k in &mut bspline.knots { *k *= tau; }
            Curve2d::BSpline(bspline)
        },
        Err(e) => {
            let avg_v = pts.iter().map(|p| p.y).sum::<f64>() / pts.len() as f64;
            Curve2d::Line(Line2d {
                origin: DVec2::new(pts[0].x, avg_v),
                direction: DVec2::new(1.0, 0.0),
            })
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cylinder functions
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the PCurve of a [`Circle3`] on a [`CylindricalSurface`].
///
/// For a circle perpendicular to the cylinder axis (axial/parallel circle),
/// returns a horizontal [`Line2d`] at v = h.  For circles in diagonal planes
/// (e.g. the Steinmetz TwoCircles case), samples 33 points on the circle,
/// maps each to cylinder UV via the analytic inverse mapping, and interpolates
/// a [`BSplineCurve2`].
pub fn circle_pcurve_on_cylinder(circle: &Circle3, cyl: &CylindricalSurface) -> Curve2d {
    let axis = cyl.axis.normalize_or_zero();
    let normal_dot = circle.normal.normalize().dot(axis).abs();
    // OCCT ProjLib_Cylinder.cxx L122-156: perp circle → Line2d with U offset and direction sign
    if (normal_dot - 1.0).abs() < TOLERANCE_MESH_LEGACY {
        let h = (circle.center - cyl.origin).dot(axis);
        // Compute U offset: angular position of circle's X-axis on cylinder
        let ux = any_perpendicular(circle.normal).normalize();
        let cyl_x = cyl.ref_dir.normalize();
        let cyl_y = axis.cross(cyl_x).normalize();
        let u = ux.dot(cyl_y).atan2(ux.dot(cyl_x));
        // Direction sign: circle normal vs cylinder axis (OCCT L145: ZCyl.Dot(aCircPos.Direction()))
        let dir_sgn = if circle.normal.dot(axis) > 0.0 { 1.0 } else { -1.0 };
        return Curve2d::Line(Line2d {
            origin: DVec2::new(u, h),
            direction: DVec2::new(dir_sgn, 0.0),
        });
    }
    // Diagonal circle: sample and map to UV.
    let u_axis = any_perpendicular(axis);
    let v_axis = axis.cross(u_axis).normalize();
    let n = 33_usize;
    let mut pts: Vec<DVec2> = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / (n - 1) as f64;
            let p3 = circle.point_at(t);
            let v = (p3 - cyl.origin).dot(axis);
            let radial = p3 - cyl.origin - axis * v;
            let mut u = radial.dot(v_axis).atan2(radial.dot(u_axis));
            if u < 0.0 {
                u += std::f64::consts::TAU;
            }
            DVec2::new(u, v)
        })
        .collect();
    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }
    match interpolate_points_2d(&pts) {
        Ok(mut bspline) => {
            // OCCT-aligned: rescale knot vector from [0, 1] to [0, TAU] to match
            // the 3D circle curve's parameter range.
            let tau = std::f64::consts::TAU;
            for k in &mut bspline.knots { *k *= tau; }
            Curve2d::BSpline(bspline)
        },
        Err(_) => {
            let a = pts[0];
            let b = *pts.last().unwrap_or(&a);
            let d = b - a;
            let dir = if d.length_squared() > TOLERANCE_VEC_SQ_MIN {
                d.normalize()
            } else {
                DVec2::X
            };
            Curve2d::Line(Line2d { origin: a, direction: dir })
        }
    }
}

/// Compute the azimuth θ of a line's origin on a cylindrical surface, mapped
/// to [0, 2π).  This is the u-coordinate the line's pcurve would have.
///
/// A line parallel to the cylinder axis at angular position θ returns θ in [0, 2π).
pub fn line_theta_on_cylinder(line: &Line3, cyl: &CylindricalSurface) -> f64 {
    let u_axis = cyl.ref_dir.normalize();
    let v_axis = cyl.axis.cross(u_axis).normalize();

    let radial = line.origin - cyl.origin;
    let radial_perp = radial - cyl.axis * radial.dot(cyl.axis.normalize());
    let mut theta = radial_perp.dot(v_axis).atan2(radial_perp.dot(u_axis));

    // Map to [0, 2π)
    if theta < 0.0 {
        theta += std::f64::consts::TAU;
    }
    theta
}

/// Compute the PCurve of a [`Line3`] on a [`CylindricalSurface`].
///
/// A line parallel to the cylinder axis at azimuth θ returns a vertical
/// [`Line2d`] at u = θ in (θ, h) space.
pub fn line_pcurve_on_cylinder(line: &Line3, cyl: &CylindricalSurface) -> Curve2d {
    let theta = line_theta_on_cylinder(line, cyl);
    let h = (line.origin - cyl.origin).dot(cyl.axis.normalize());

    Curve2d::Line(Line2d {
        origin: DVec2::new(theta, h),
        direction: DVec2::new(0.0, 1.0),
    })
}

/// Compute the PCurve of an [`Ellipse3`] on a [`CylindricalSurface`].
///
/// Samples the ellipse at 33 evenly-spaced parameter values over [0, 2π],
/// projects each 3D point onto the cylinder's (u, v) domain using the analytic
/// inverse mapping, and interpolates a [`BSplineCurve2`].
///
/// Unlike [`fallback_pcurve_by_projection`], this avoids iterative closest-point
/// projection and is exact for points on the cylinder surface. The ellipse's
/// own parameterization gives well-distributed sample points.
pub fn ellipse_pcurve_on_cylinder(ellipse: &Ellipse3, cyl: &CylindricalSurface) -> Curve2d {
    let axis = cyl.axis.normalize_or_zero();
    let u_axis = cyl.ref_dir.normalize();
    let v_axis = axis.cross(u_axis).normalize();

    let n = 33_usize;
    let mut pts: Vec<DVec2> = (0..n)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / (n - 1) as f64;
            let p3 = ellipse.point_at(t);
            let v = (p3 - cyl.origin).dot(axis);
            let radial = p3 - cyl.origin - axis * v;
            let mut u = radial.dot(v_axis).atan2(radial.dot(u_axis));
            // Map to [0, 2π] to match cylinder UV convention (same as
            // closest_point_on_surface and cone_uv_from_point).
            if u < 0.0 {
                u += std::f64::consts::TAU;
            }
            DVec2::new(u, v)
        })
        .collect();

    // Unwrap seam discontinuities: samples near the 0/2π seam may have a
    // large jump if the ellipse crosses it.
    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }

    match interpolate_points_2d(&pts) {
        Ok(mut bspline) => {
            // OCCT-aligned: rescale knot vector from [0, 1] to [0, TAU] to match
            // the 3D ellipse curve's parameter range.
            let tau = std::f64::consts::TAU;
            for k in &mut bspline.knots { *k *= tau; }
            Curve2d::BSpline(bspline)
        },
        Err(_) => {
            let a = pts[0];
            let b = *pts.last().unwrap_or(&a);
            let d = b - a;
            let dir = if d.length_squared() > TOLERANCE_VEC_SQ_MIN {
                d.normalize()
            } else {
                DVec2::X
            };
            Curve2d::Line(Line2d { origin: a, direction: dir })
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cone functions
// ─────────────────────────────────────────────────────────────────────────────

fn cone_uv_from_point(point: DVec3, cone: &ConicalSurface) -> DVec2 {
    let axis = cone.axis_dir();
    let u_axis = any_perpendicular(axis);
    let v_axis = axis.cross(u_axis).normalize();
    let local = point - cone.apex;
    let axial = local.dot(axis);
    let radial = local - axis * axial;
    let mut u = radial.dot(v_axis).atan2(radial.dot(u_axis));
    if u < 0.0 {
        u += std::f64::consts::TAU;
    }
    DVec2::new(u, cone.slant_from_axial(axial))
}

fn sampled_curve_pcurve_on_cone(
    curve: &rcad_kernel::geom::Curve3,
    t_range: &[f64; 2],
    cone: &ConicalSurface,
) -> Curve2d {
    let n = 33_usize;
    let mut pts: Vec<DVec2> = (0..n)
        .map(|i| {
            let t = t_range[0] + (t_range[1] - t_range[0]) * i as f64 / (n - 1) as f64;
            let p3 = curve.point_at(t);
            cone_uv_from_point(p3, cone)
        })
        .collect();

    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }

    let mut bspline = interpolate_points_2d(&pts).expect("cone curve samples should not be degenerate");
    // OCCT-aligned: rescale knot vector from [0, 1] to [t_range[0], t_range[1]] to match
    // the 3D curve's parameter range.
    let ts = t_range[0];
    let te = t_range[1];
    let span = te - ts;
    if span > 0.0 {
        for k in &mut bspline.knots { *k = ts + (*k) * span; }
    }
    Curve2d::BSpline(bspline)
}

pub fn circle_pcurve_on_cone(circle: &Circle3, cone: &ConicalSurface) -> Curve2d {
    let axis = cone.axis_dir();
    let normal_dot = circle.normal.normalize().dot(axis).abs();
    if (normal_dot - 1.0).abs() < TOLERANCE_MESH_LEGACY {
        let slant = cone.slant_from_axial((circle.center - cone.apex).dot(axis));
        return Curve2d::Line(Line2d {
            origin: DVec2::new(0.0, slant),
            direction: DVec2::new(1.0, 0.0),
        });
    }
    sampled_curve_pcurve_on_cone(&rcad_kernel::geom::Curve3::Circle(*circle), &[0.0, std::f64::consts::TAU], cone)
}

pub fn line_pcurve_on_cone(line: &Line3, cone: &ConicalSurface) -> Curve2d {
    let uv0 = cone_uv_from_point(line.origin, cone);
    let uv1 = cone_uv_from_point(line.origin + line.direction, cone);
    let du = (uv1.x - uv0.x).abs().min((uv1.x - uv0.x + std::f64::consts::TAU).abs());
    if du < TOLERANCE_MESH_LEGACY {
        let dir_v = if uv1.y >= uv0.y { 1.0 } else { -1.0 };
        return Curve2d::Line(Line2d {
            origin: uv0,
            direction: DVec2::new(0.0, dir_v),
        });
    }
    sampled_curve_pcurve_on_cone(&rcad_kernel::geom::Curve3::Line(*line), &[-10.0, 10.0], cone)
}

pub fn ellipse_pcurve_on_cone(ellipse: &Ellipse3, cone: &ConicalSurface) -> Curve2d {
    sampled_curve_pcurve_on_cone(&rcad_kernel::geom::Curve3::Ellipse(*ellipse), &[0.0, std::f64::consts::TAU], cone)
}

pub fn sampled_pcurve_on_cone(
    curve: &rcad_kernel::geom::Curve3,
    t_range: &[f64; 2],
    cone: &ConicalSurface,
) -> Curve2d {
    sampled_curve_pcurve_on_cone(curve, t_range, cone)
}

// ─────────────────────────────────────────────────────────────────────────────
// Numeric fallback functions
// ─────────────────────────────────────────────────────────────────────────────

/// Derive a PCurve by sampling `curve` at 33 evenly-spaced parameter values
/// over `t_range` and projecting each 3D point onto `surface`.
///
/// Intended as a fallback for curve/surface combinations that do not have an
/// analytic form.  Returns a [`BSplineCurve2`] interpolated through the
/// projected (u, v) points, wrapped in [`TrimmedCurve2`] to preserve the
/// mapping between the 3D curve's parameter range and the BSpline's native
/// `[0, 1]` parameterization.
///
/// When projection collapses (e.g. very short edges), fall back to a UV line.
pub fn fallback_pcurve_by_projection(
    curve: &rcad_kernel::geom::Curve3,
    t_range: &[f64; 2],
    surface: &Surface3,
) -> Curve2d {
    let n = 33_usize;
    let mut pts: Vec<DVec2> = (0..n)
        .map(|i| {
            let t = t_range[0] + (t_range[1] - t_range[0]) * i as f64 / (n - 1) as f64;
            let p3 = curve.point_at(t);
            match surface {
                Surface3::Sphere(sph) => sph.world_to_uv(p3),
                Surface3::Cone(cone) => cone.world_to_uv(p3),
                Surface3::Torus(torus) => torus.world_to_uv(p3),
                _ => {
                    let proj = closest_point_on_surface(surface, p3, 16);
                    DVec2::new(proj.params.0, proj.params.1)
                }
            }
        })
        .collect();

    // Unwrap seam discontinuities: u values are in [0, 2π] after the mapping,
    // but consecutive samples may still jump by ~2π when the ellipse crosses
    // the seam. Make the u sequence monotone for a clean BSpline.
    // Make the u sequence monotone so the interpolated BSpline has no kinks.
    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            // Jumped from near -π back up to near +π: pull remaining down.
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            // Jumped from near +π down to near -π: push remaining up.
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }

    match interpolate_points_2d(&pts) {
        Ok(bspline) => {
            // ✅ OCCT对齐: Geom2d_TrimmedCurve — wrap in TrimmedCurve2 to
            // preserve the mapping between 3D curve parameter range [t0, t1]
            // and BSpline's native [0, 1] parameterization.
            let tc = rcad_kernel::geom::TrimmedCurve2 {
                curve: Box::new(Curve2d::BSpline(bspline)),
                t_min: t_range[0],
                t_max: t_range[1],
            };
            Curve2d::Trimmed(tc)
        }
        Err(_) => {
            let a = pts[0];
            let b = *pts.last().unwrap_or(&a);
            let d = b - a;
            let dir = if d.length_squared() > TOLERANCE_VEC_SQ_MIN {
                d.normalize()
            } else {
                DVec2::X
            };
            Curve2d::Line(Line2d { origin: a, direction: dir })
        }
    }
}

/// Project a 3D polyline onto `surface` and interpolate a [`BSplineCurve2`].
///
/// Returns `None` if the polyline has fewer than 2 points or all projected
/// points are coincident.
pub fn polyline_pcurve_by_projection(polyline: &[DVec3], surface: &Surface3) -> Option<Curve2d> {
    if polyline.len() < 2 {
        return None;
    }

    let mut pts: Vec<DVec2> = polyline
        .iter()
        .map(|&p3| match surface {
            Surface3::Sphere(sph) => sph.world_to_uv(p3),
            _ => {
                let proj = closest_point_on_surface(surface, p3, 16);
                DVec2::new(proj.params.0, proj.params.1)
            }
        })
        .collect();

    // Phase 1: Standard seam unwrap using π-threshold.
    if std::env::var("RCAD_DEBUG_PC").is_ok() && matches!(surface, Surface3::Cylinder(_)) {
        eprintln!("[PC_DEBUG] Phase1 start: pts[0].x={:.6} pts[last].x={:.6} n={}",
            pts[0].x, pts[pts.len()-1].x, pts.len());
    }
    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }

    if std::env::var("RCAD_DEBUG_PC").is_ok() && matches!(surface, Surface3::Cylinder(_)) {
        eprintln!("[PC_DEBUG] Phase1 end: pts[0].x={:.6} pts[last].x={:.6}", pts[0].x, pts[pts.len()-1].x);
        if pts.len() >= 5 {
            eprintln!("[PC_DEBUG] Phase1 sample: [{:.6}, {:.6}, {:.6}, {:.6}, {:.6}]",
                pts[0].x, pts[1].x, pts[2].x, pts[3].x, pts[4].x);
        }
    }

    // Phase 2: Detect V-shape folds. When the polyline has two disconnected curve
    // segments concatenated (e.g. PerpendicularOffsetCurves with a near-π gap),
    // the standard π-threshold unwrap may fold the sequence — u goes up then back
    // instead of continuing in one direction.
    //
    // Detection: 1 sign change in consecutive du values (e.g. all-negative →
    // all-positive), with span < 2π and near-zero net delta.
    //
    // NOTE: On Cylinder surfaces, this V-fold is geometrically correct — the
    // analytic atan2 formula gives identical values to Newton projection, and the
    // intersection curve genuinely wraps from u=0 to u=π and back. The fix for
    // offset-cylinder boolean failures (ZE7-9, ZF1-4) is in the Pave-Filler, not
    // here. This detection is diagnostic-only.
    if pts.len() >= 3 && matches!(surface, Surface3::Cylinder(_)) {
        let u_min = pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let u_max = pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let abs_span = u_max - u_min;
        let net_delta = pts[pts.len() - 1].x - pts[0].x;
        if abs_span < std::f64::consts::TAU * 0.9 && net_delta.abs() < abs_span * 0.5 {
            // Count sign changes in consecutive du values to detect V-shape fold.
            // A V-shape fold has exactly one sign change (e.g. all-negative → all-positive).
            // A boundary jump between valid segments has two sign changes (e.g.
            // positive → single-negative-jump → positive).
            let mut sign_changes = 0;
            let mut fold_idx = 0;
            let mut prev_sign = 0i8;
            for i in 1..pts.len() {
                let du = pts[i].x - pts[i - 1].x;
                let sign = if du.abs() > 1e-12 {
                    if du > 0.0 { 1 } else { -1 }
                } else { 0 };
                if sign != 0 && prev_sign != 0 && sign != prev_sign {
                    sign_changes += 1;
                    if sign_changes == 1 {
                        fold_idx = i;
                    }
                }
                if sign != 0 {
                    prev_sign = sign;
                }
            }
            eprintln!("[PCURVE_SIGNS] n_pts={} abs_span={:.6} net_delta={:.6} sign_changes={} fold_idx={}",
                pts.len(), abs_span, net_delta, sign_changes, fold_idx);
            if sign_changes == 1 && fold_idx >= 2 && fold_idx + 2 <= pts.len() {
                eprintln!("[PCURVE_FOLD] V-fold: span={:.6} net={:.6} fold_idx={} n={} u_first={:.4} u_fold={:.4} u_last={:.4}",
                    abs_span, net_delta, fold_idx, pts.len(),
                    pts[0].x, pts[fold_idx].x, pts[pts.len()-1].x);
                if let Surface3::Cylinder(cyl) = surface {
                    // Analytic u recomputation is NOT effective here — the V-fold
                    // is geometrically correct on the cylinder UV (the intersection
                    // curve wraps around the front of the cylinder and back). The
                    // Newton projection already gives the correct analytic u values
                    // (identical to atan2). The V-fold detection is diagnostic only.
                    //
                    // The actual fix for offset cylinder failures (ZE7-9, ZF1-4) is
                    // in the Pave-Filler: the PerpendicularOffsetCurves handler
                    // concatenates two ~π segments with a near-π gap, and Phase 1
                    // unwrap can choose the wrong branch. The BSpline through these
                    // folded points produces a self-intersecting PCurve, which causes
                    // the boolean builder to produce wrong geometry.
                    let _ = cyl;
                    eprintln!("[PCURVE_FOLD] V-fold on Cylinder — diagnostic only, no fix applied");
                }
                // Note: for Cone/Torus surfaces we do not apply a fold fix
                // here.  Empirical testing shows the V-shape fold only occurs
                // on cylinder surfaces from Newton's 2π-periodic branch issue.
            }
        }
    }

    interpolate_points_2d(&pts).ok().map(Curve2d::BSpline)
}

/// Check whether `pcurve` lies entirely within the given UV bounds.
///
/// ✅ OCCT对齐: CheckPCurve (IntTools_FaceFace.cxx L2924-2999)
///
/// Samples each C0-interval at NPoints (OCCT uses 23 per interval) and
/// verifies every sample lies within the UV bounds (plus a relative
/// tolerance of 1% of span, minimum `TOLERANCE_ABS`).
///
/// For periodic surfaces, shifts the UV bounds by whole periods so the
/// midpoint of the pcurve falls within the shifted domain.  This matches
/// OCCT's approach of evaluating the midpoint first and shifting accordingly.
///
/// `t_range` is the effective parameter range to sample — for most curves
/// this is the 3D curve's [`IntersectionCurve::t_range`]; for BSpline/Bezier
/// pcurves it should be `[0.0, 1.0]`.
///
/// Returns `true` when all sampled points are within bounds (or the curve
/// is degenerate / has no finite range).
pub fn check_pcurve_in_face(
    pcurve: &Curve2d,
    t_range: [f64; 2],
    uv_bounds: [f64; 4],
    u_period: Option<f64>,
    v_period: Option<f64>,
) -> bool {
    const N_POINTS: usize = 23;

    let [umin, umax, vmin, vmax] = uv_bounds;
    let tol_u = ((umax - umin) * 0.01).max(TOLERANCE_ABS);
    let tol_v = ((vmax - vmin) * 0.01).max(TOLERANCE_ABS);

    let [t0, t1] = t_range;
    if !t0.is_finite() || !t1.is_finite() || (t1 - t0).abs() < TOLERANCE_LEN_MIN {
        return true; // degenerate range — skip
    }

    // Periodic shift: shift UV bounds so the midpoint pcurve parameter
    // falls in the shifted domain (matching OCCT's approach).
    let mid_t = 0.5 * (t0 + t1);
    let mid_uv = match pcurve {
        Curve2d::Trimmed(tc) => {
            // For a TrimmedCurve2, the effective range is [t_min, t_max],
            // but point_at(t) for t in that range maps correctly.  Compute
            // the midpoint in the inner curve's native range by mapping.
            let mt = 0.5 * (tc.t_min + tc.t_max);
            tc.curve.as_ref().point_at(mt)
        }
        other => other.point_at(mid_t),
    };

    let u_shift = u_period
        .map(|per| {
            let raw = mid_uv.x - umin;
            (raw / per).floor() * per
        })
        .unwrap_or(0.0);
    let v_shift = v_period
        .map(|per| {
            let raw = mid_uv.y - vmin;
            (raw / per).floor() * per
        })
        .unwrap_or(0.0);

    // Sample N_POINTS evenly over the parameter range
    for i in 0..N_POINTS {
        let t = t0 + (t1 - t0) * i as f64 / (N_POINTS - 1) as f64;
        let uv = pcurve.point_at(t);
        let u = uv.x - u_shift;
        let v = uv.y - v_shift;
        if umin - u > tol_u || u - umax > tol_u || vmin - v > tol_v || v - vmax > tol_v {
            return false;
        }
    }
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// IsCurveValid — pcurve self-intersection check
// ─────────────────────────────────────────────────────────────────────────────

/// Check whether a 2D curve is free of self-intersections.
///
/// ✅ OCCT对齐: IsCurveValid (IntTools_FaceFace.cxx L2252-2289)
///
/// For analytic curves (Line, Circle, Ellipse, etc.) trivially returns true.
/// For BSpline/Bezier curves, samples the curve and checks for non-adjacent
/// segment intersections.
///
/// Returns `false` if the curve is self-intersecting or null.
pub fn is_curve_valid_2d(curve: &Curve2d) -> bool {
    match curve {
        Curve2d::Trimmed(tc) => is_curve_valid_2d(tc.curve.as_ref()),
        // Analytic curves never self-intersect
        Curve2d::Line(_)
        | Curve2d::Circle(_)
        | Curve2d::Ellipse(_)
        | Curve2d::CircleInvolute(_)
        | Curve2d::ArchimedeanSpiral(_)
        | Curve2d::LogarithmicSpiral(_)
        | Curve2d::SineWave(_)
        | Curve2d::Parabola(_)
        | Curve2d::Hyperbola(_) => true,
        // BSpline/Bezier: check polyline self-intersection
        Curve2d::BSpline(_) | Curve2d::Bezier(_) => {
            let pts = (0..100)
                .map(|i| {
                    let t = i as f64 / 99.0; // [0, 1]
                    curve.point_at(t)
                })
                .collect::<Vec<_>>();
            !has_polyline_self_intersection_2d(&pts)
        }
    }
}

/// Check if a polyline has non-adjacent segments that intersect.
fn has_polyline_self_intersection_2d(points: &[DVec2]) -> bool {
    if points.len() < 4 {
        return false;
    }
    for i in 0..points.len() - 1 {
        for j in (i + 2)..points.len() - 1 {
            if segments_intersect_2d_open(points[i], points[i + 1], points[j], points[j + 1]) {
                return true;
            }
        }
    }
    false
}

/// Check if two 2D open segments intersect in their interiors.
/// Returns false for parallel, collinear, or endpoint-only intersections.
fn segments_intersect_2d_open(p1: DVec2, p2: DVec2, p3: DVec2, p4: DVec2) -> bool {
    let d1 = p2 - p1;
    let d2 = p4 - p3;
    let cross = d1.perp_dot(d2);
    if cross.abs() < TOLERANCE_LEN_MIN {
        return false; // parallel
    }
    let dp = p3 - p1;
    let t = dp.perp_dot(d2) / cross;
    let u = dp.perp_dot(d1) / cross;
    // Open interval (0, 1) — endpoints are shared vertices
    t > TOLERANCE_ABS && t < 1.0 - TOLERANCE_ABS && u > TOLERANCE_ABS && u < 1.0 - TOLERANCE_ABS
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeTolReached3d — deviation between 3D curve and pcurve
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the maximum deviation between a 3D curve and its pcurve on a surface.
///
/// ✅ OCCT对齐: IntTools_Tools::ComputeTolerance / FindMaxDistance
/// (IntTools_FaceFace.cxx L603-681, L2813-2918)
///
/// Uses golden-section search to find `t` in `[t0, t1]` that maximises
/// `||C3D(t) - surface(pcurve(t))||`.
///
/// Returns the maximum deviation in model units.
pub fn compute_max_deviation_3d_to_pcurve(
    curve_3d: &rcad_kernel::geom::Curve3,
    pcurve: &Curve2d,
    surface: &Surface3,
    t_range: [f64; 2],
) -> f64 {
    let [t0, t1] = t_range;
    if !t0.is_finite() || !t1.is_finite() || (t1 - t0).abs() < TOLERANCE_LEN_MIN {
        return 0.0;
    }
    let f = |t: f64| {
        let p3 = curve_3d.point_at(t);
        let uv = pcurve.point_at(t);
        let p_surf = surface.point_at(uv.x, uv.y);
        (p3 - p_surf).length()
    };
    crate::golden_section_max(f, t0, t1, TOLERANCE_PARAM_LEGACY)
}

/// Compute the maximum deviation between a 3D curve and a surface
/// (when no pcurve is available).
///
/// ✅ OCCT对齐: FindMaxDistance (IntTools_FaceFace.cxx L2813-2847)
///
/// Projects equally-spaced samples onto the surface and uses golden-section
/// search per segment to find the maximum distance.
pub fn compute_max_deviation_from_surface(
    curve_3d: &rcad_kernel::geom::Curve3,
    surface: &Surface3,
    t_range: [f64; 2],
) -> f64 {
    let [t0, t1] = t_range;
    if !t0.is_finite() || !t1.is_finite() || (t1 - t0).abs() < TOLERANCE_LEN_MIN {
        return 0.0;
    }
    // Divide into 11 segments (OCCT uses aNbS = 11)
    let n_seg = 11_usize;
    let dt = (t1 - t0) / n_seg as f64;
    let an_eps = 1e-4 * dt;
    let mut max_d = 0.0;
    for seg in 0..n_seg {
        let seg_start = t0 + seg as f64 * dt;
        let seg_end = (seg_start + dt).min(t1);
        let f = |t: f64| {
            let p3 = curve_3d.point_at(t);
            let proj = closest_point_on_surface(surface, p3, 16);
            (p3 - proj.point).length()
        };
        let d = crate::golden_section_max(f, seg_start, seg_end, an_eps);
        if d > max_d {
            max_d = d;
        }
    }
    max_d
}

/// Compute the tolerance and tangential tolerance for an intersection curve
/// by evaluating the deviation between its 3D curve and pcurves on both surfaces.
///
/// ✅ OCCT对齐: ComputeTolReached3d (IntTools_FaceFace.cxx L603-681)
///
/// `current_tol` is the starting tolerance (e.g. from the intersection algorithm).
/// Returns `(updated_tolerance, tangential_tolerance)`.
pub fn compute_intersection_curve_tolerance(
    curve_3d: &rcad_kernel::geom::Curve3,
    pcurve_on_a: Option<&Curve2d>,
    pcurve_on_b: Option<&Curve2d>,
    surface_a: &Surface3,
    surface_b: &Surface3,
    t_range: [f64; 2],
    face_tol_a: f64,
    face_tol_b: f64,
    current_tol: f64,
) -> (f64, f64) {
    let mut tol = current_tol;
    let [t0, t1] = t_range;
    // PCurve on surface A
    if let Some(pca) = pcurve_on_a {
        let d = compute_max_deviation_3d_to_pcurve(curve_3d, pca, surface_a, [t0, t1]);
        if d > tol {
            tol = d;
        }
    } else {
        let d = compute_max_deviation_from_surface(curve_3d, surface_a, [t0, t1]);
        if d > tol {
            tol = d;
        }
    }
    // PCurve on surface B
    if let Some(pcb) = pcurve_on_b {
        let d = compute_max_deviation_3d_to_pcurve(curve_3d, pcb, surface_b, [t0, t1]);
        if d > tol {
            tol = d;
        }
    } else {
        let d = compute_max_deviation_from_surface(curve_3d, surface_b, [t0, t1]);
        if d > tol {
            tol = d;
        }
    }
    // Tangential tolerance: at least the max face tolerance
    let tang_tol = face_tol_a.max(face_tol_b);
    (tol, tang_tol)
}

// ─────────────────────────────────────────────────────────────────────────────
// PrepareLines3D — closed‑curve splitting + redundant‑line filtering
// ─────────────────────────────────────────────────────────────────────────────

/// Split a closed intersection curve (full circle / ellipse) into two halves
/// at the u-parameter midpoint.
///
/// ✅ OCCT对齐: IntTools_Tools::SplitCurve (用于 PrepareLines3D)
///
/// When a closed curve has the full [0, 2π] parametric range it cannot be
/// properly trimmed as a single segment — OCCT splits it into complementary
/// arcs.  Returns `None` for curves that are not closed or not full-range.
fn split_closed_curve(
    curve_3d: &rcad_kernel::geom::Curve3,
    t_range: &[f64; 2],
) -> Option<[[f64; 2]; 2]> {
    let [t0, t1] = *t_range;
    let is_full_circle = match curve_3d {
        rcad_kernel::geom::Curve3::Circle(_) | rcad_kernel::geom::Curve3::Ellipse(_) => {
            (t1 - t0 - std::f64::consts::TAU).abs() < TOLERANCE_ANG
        }
        _ => false,
    };
    if !is_full_circle {
        return None;
    }
    // Split at mid-point of the range
    let tm = 0.5 * (t0 + t1);
    Some([[t0, tm], [tm, t1]])
}

/// Post-process intersection curves: split closed curves and reject redundant
/// lines.
///
/// ✅ OCCT对齐: PrepareLines3D (IntTools_FaceFace.cxx L1898-1979)
///
/// 1. Splits closed 3D curves (circles/ellipses with full [0, 2π] range) at
///    the parametric midpoint so they can be trimmed properly.
/// 2. (Future) Plane/Cone 4-line redundant-line rejection.
///
/// Operates on the `curves` vector in place.
pub fn prepare_lines_3d(curves: &mut Vec<crate::bopds::ds::IntersectionCurve>) {
    let mut new_curves: Vec<crate::bopds::ds::IntersectionCurve> = Vec::new();

    for ic in curves.drain(..) {
        // 1. Split closed curves
        let splits = split_closed_curve(&ic.curve, &ic.t_range);
        if let Some([r0, r1]) = splits {
            let mut c0 = ic.clone();
            c0.t_range = r0;
            let mut c1 = ic;
            c1.t_range = r1;
            new_curves.push(c0);
            new_curves.push(c1);
        } else {
            new_curves.push(ic);
        }
    }

    // 2. (Plane/Cone 4-line rejection — reserved for future alignment)

    *curves = new_curves;
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// Adaptive polyline refinement
// ─────────────────────────────────────────────────────────────────────────────

/// Refine a polyline by subdividing segments where chord error exceeds tolerance.
///
/// `samples` is a slice of `(parameter, point_3d)` pairs representing the curve
/// at discrete parameter values.  For each consecutive pair, evaluates the curve
/// at the parametric midpoint using `eval_fn`.  If the midpoint deviates from the
/// chord by more than `chord_tol`, the segment is recursively subdivided up to
/// `max_depth` levels.
///
/// Returns a new `Vec<(f64, DVec3)>` with additional inserted samples.
pub fn refine_polyline<F>(
    samples: &[(f64, DVec3)],
    eval_fn: F,
    chord_tol: f64,
    max_depth: usize,
) -> Vec<(f64, DVec3)>
where
    F: Fn(f64) -> Option<DVec3>,
{
    if samples.len() < 2 {
        return samples.to_vec();
    }

    let chord_tol_sq = chord_tol * chord_tol;

    fn subdivide<F>(
        p0: DVec3,
        p1: DVec3,
        u0: f64,
        u1: f64,
        eval_fn: &F,
        chord_tol_sq: f64,
        depth: usize,
        max_depth: usize,
        out: &mut Vec<(f64, DVec3)>,
    ) where
        F: Fn(f64) -> Option<DVec3>,
    {
        if depth >= max_depth {
            return;
        }

        let u_mid = (u0 + u1) * 0.5;
        if let Some(p_mid) = eval_fn(u_mid) {
            let chord = p1 - p0;
            let chord_len_sq = chord.length_squared();
            if chord_len_sq > 0.0 {
                let t = ((p_mid - p0).dot(chord) / chord_len_sq).clamp(0.0, 1.0);
                let chord_pt = p0 + t * chord;
                if (p_mid - chord_pt).length_squared() > chord_tol_sq {
                    // Subdivide left and right
                    subdivide(p0, p_mid, u0, u_mid, eval_fn, chord_tol_sq, depth + 1, max_depth, out);
                    out.push((u_mid, p_mid));
                    subdivide(p_mid, p1, u_mid, u1, eval_fn, chord_tol_sq, depth + 1, max_depth, out);
                }
            }
        }
    }

    let mut result = Vec::with_capacity(samples.len() * 2);
    result.push(samples[0]);

    for i in 0..samples.len() - 1 {
        let (u0, p0) = samples[i];
        let (u1, p1) = samples[i + 1];

        subdivide(p0, p1, u0, u1, &eval_fn, chord_tol_sq, 0, max_depth, &mut result);

        // Append the end point if it was not already inserted by subdivision
        let last_u = result.last().unwrap().0;
        if (last_u - u1).abs() > 1e-14 {
            result.push((u1, p1));
        }
    }

    // Trim trailing near-duplicates (common for closed curves)
    while result.len() >= 3 {
        let n = result.len();
        let d = (result[n - 1].1 - result[0].1).length_squared();
        if d < TOLERANCE_VEC_SQ_MIN {
            result.pop();
        } else {
            break;
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Curve2dEval, SphericalSurface};
    use std::f64::consts::PI;

    /// A circle whose normal is Z lying in the XY plane (z = 0) projects to a
    /// Circle2d in the plane's (u, v) space.
    #[test]
    fn circle_on_plane_is_circle() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let circle = Circle3::new(DVec3::new(1.0, DVec3::Z, 3.0);

        let pcurve = circle_pcurve_on_plane(&circle, &plane);

        match pcurve {
            Curve2d::Circle(c) => {
                assert!((c.radius - 3.0).abs() < TOLERANCE_COORD_SUB, "radius={}", c.radius);
            }
            other => panic!("expected Circle2d, got {other:?}"),
        }
    }

    /// A circle at z = 1 on a sphere of radius 2 (axis = Z, center = origin)
    /// should produce φ = acos(0.5) = π/3.
    #[test]
    fn circle_on_sphere_is_latitude() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Z),
        };
        let circle = Circle3::new(DVec3::new(0.0, DVec3::Z, (3.0_f64).sqrt());

        let pcurve = circle_pcurve_on_sphere(&circle, &sphere);

        match pcurve {
            Curve2d::Line(l) => {
                let expected_phi = (0.5_f64).acos(); // π/3
                assert!(
                    (l.origin.y - expected_phi).abs() < TOLERANCE_COORD_SUB,
                    "phi={}, expected {expected_phi}",
                    l.origin.y
                );
                // Origin x starts at -π so the line spans [-π, +π] over [0, 2π] sampling.
                assert!(
                    (l.origin.x + PI).abs() < TOLERANCE_COORD_SUB,
                    "expected origin.x = -π, got {}",
                    l.origin.x
                );
                // Direction must be horizontal (constant colatitude).
                assert!((l.direction.x - 1.0).abs() < TOLERANCE_COORD_SUB);
                assert!(l.direction.y.abs() < TOLERANCE_COORD_SUB);
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    /// A circle at height h = 3 on a cylinder should produce a horizontal
    /// line at v = 3.
    #[test]
    fn circle_on_cylinder_is_h_line() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        };
        let circle = Circle3::new(DVec3::new(0.0, DVec3::Z, 1.0);

        let pcurve = circle_pcurve_on_cylinder(&circle, &cyl);

        match pcurve {
            Curve2d::Line(l) => {
                assert!(
                    (l.origin.y - 3.0).abs() < TOLERANCE_COORD_SUB,
                    "h={}, expected 3.0",
                    l.origin.y
                );
                assert!((l.direction.x - 1.0).abs() < TOLERANCE_COORD_SUB);
                assert!(l.direction.y.abs() < TOLERANCE_COORD_SUB);
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    #[test]
    fn circle_on_cone_is_h_line() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: (1.0_f64 / 2.0).atan(),
        };
        let h = 3.0;
        let slant = cone.slant_from_axial(h);
        let circle = Circle3::new(DVec3::new(0.0, DVec3::Z, h * cone.half_angle_rad.tan());

        let pcurve = circle_pcurve_on_cone(&circle, &cone);
        match pcurve {
            Curve2d::Line(l) => {
                assert!((l.origin.y - slant).abs() < TOLERANCE_COORD_SUB);
                assert!((l.direction.x - 1.0).abs() < TOLERANCE_COORD_SUB);
                assert!(l.direction.y.abs() < TOLERANCE_COORD_SUB);
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    #[test]
    fn line_on_cone_is_v_line() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: (1.0_f64 / 2.0).atan(),
        };
        let line = Line3 {
            origin: DVec3::new(1.0, 0.0, 2.0),
            direction: DVec3::new(0.5, 0.0, 1.0).normalize(),
        };

        let pcurve = line_pcurve_on_cone(&line, &cone);
        match pcurve {
            Curve2d::Line(l) => {
                // The origin's u coordinate depends on the arbitrary perpendicular chosen,
                // so we only verify that the line is a v-line (direction purely in v)
                // by checking that the x direction is zero.
                assert!(l.direction.x.abs() < TOLERANCE_COORD_SUB, "v-line should have zero u direction");
                assert!((l.direction.y - 1.0).abs() < TOLERANCE_COORD_SUB, "v-line should have unit v direction");
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    /// The fallback projection of any curve on a sphere should produce a
    /// BSplineCurve2.
    #[test]
    fn fallback_projection_produces_bspline() {
        use rcad_kernel::geom::Curve3;

        let sphere_surface = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        // A circle on the sphere at the equator (z = 0, r = 2).
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0,
        );
        let curve3 = Curve3::Circle(circle);
        let t_range = [0.0_f64, PI]; // half circle

        let pcurve = fallback_pcurve_by_projection(&curve3, &t_range, &sphere_surface);

        match pcurve {
            Curve2d::BSpline(ref b) => {
                // Should have at least some control points.
                assert!(!b.control_points.is_empty());
                // Evaluate endpoints to make sure the BSpline is usable.
                let p0 = pcurve.point_at(0.0);
                let p1 = pcurve.point_at(1.0);
                // Both must be finite.
                assert!(p0.x.is_finite() && p0.y.is_finite());
                assert!(p1.x.is_finite() && p1.y.is_finite());
            }
            other => panic!("expected BSpline2, got {other:?}"),
        }
    }

    /// `circle_pcurve_on_sphere` is valid for any circle that lies on the sphere,
    /// even when the circle's normal is NOT parallel to the sphere axis.
    /// Here the intersection of two spheres whose centres are separated along X
    /// gives a circle whose normal is X, but we still get a correct latitude line.
    #[test]
    fn circle_on_sphere_non_axis_normal() {
        // Two unit spheres: sph1 at origin (axis=Z), sph2 at (1,0,0).
        // Their intersection circle: d=1, r1=r2=1.
        //   h = (1 + 1 - 1)/(2) = 0.5   (distance from sph1 center to radical plane)
        //   r_circ = sqrt(1 - 0.25) = sqrt(0.75)
        //   circle center = (0.5, 0, 0)
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 1.0);
        let circle = Circle3::new(DVec3::new(0.5, DVec3::X, (0.75_f64).sqrt());

        let pcurve = circle_pcurve_on_sphere(&circle, &sphere);

        // along_axis = (0.5, 0, 0) · (0, 0, 1) = 0  →  phi = acos(0) = π/2
        match pcurve {
            Curve2d::Line(l) => {
                let expected_phi = std::f64::consts::PI / 2.0;
                assert!(
                    (l.origin.y - expected_phi).abs() < TOLERANCE_COORD_SUB,
                    "phi={}, expected π/2",
                    l.origin.y
                );
                assert!((l.direction.x - 1.0).abs() < TOLERANCE_COORD_SUB);
                assert!(l.direction.y.abs() < TOLERANCE_COORD_SUB);
                // origin.x = longitude of circle.point_at(0) — just check it's finite
                assert!(l.origin.x.is_finite());
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    /// Verify that `circle_pcurve_on_sphere` and `fallback_pcurve_by_projection`
    /// agree at the equatorial circle (both should give v ≈ π/2).
    #[test]
    fn analytic_sphere_pcurve_matches_fallback() {
        use rcad_kernel::geom::{Curve2dEval, Curve3};

        let sphere_surf = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 2.0);
        // Equatorial circle at z=0, r=2
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0 );

        let analytic = circle_pcurve_on_sphere(&circle, &sphere);
        let fallback = fallback_pcurve_by_projection(
            &Curve3::Circle(circle),
            &[0.0, std::f64::consts::TAU],
            &sphere_surf,
        );

        // Both should yield v ≈ π/2 everywhere (equator = colatitude π/2)
        for i in 0..8 {
            let t = i as f64 / 8.0;
            let pa = analytic.point_at(t);
            let pf = fallback.point_at(t);
            assert!(
                (pa.y - pf.y).abs() < 0.02,
                "t={t}: analytic.v={:.4} fallback.v={:.4}",
                pa.y,
                pf.y
            );
        }
    }
}
