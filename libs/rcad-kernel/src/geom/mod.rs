use glam::{DVec2, DVec3};

pub type Point3 = DVec3;
pub type Vec3 = DVec3;
pub type Point2 = DVec2;
pub type Vec2 = DVec2;

// Rule 5 file-size split: the payload types and the evaluation traits moved to
// child modules; every public path is preserved by the re-exports below
// (`rcad_kernel::geom::X` resolves exactly as before the split).
mod eval_traits;
mod types_curve2d;
mod types_curve3;
mod types_surface;

pub use eval_traits::*;
pub use types_curve2d::*;
pub use types_curve3::*;
pub use types_surface::*;

/// OCCT-aligned: returns a vector perpendicular to `v` using the gp_Ax2
/// reference direction convention (project X onto plane; fallback to Z).
/// Stable for any non-zero input.
pub fn any_perpendicular(v: DVec3) -> DVec3 {
    // OCCT gp_Ax2 reference direction selection:
    // Use X (1,0,0) as default; if |v·X| ≥ 1-1e-12 (v parallel to X), use Z (0,0,1).
    let ref_dir = if v.x.abs() > 1.0 - 1e-12 {
        DVec3::Z
    } else {
        DVec3::X
    };
    // Project reference onto plane perpendicular to v: ref - v*(v·ref)
    let perp = ref_dir - v * ref_dir.dot(v);
    perp.normalize_or_zero()
}

/// OCCT-aligned: 90° counter-clockwise rotation in 2D (gp_Dir2d::Rotated).
/// Equivalent to cross product with (0, 0, 1) in 2D homogeneous form.
pub fn turn_2d(v: DVec2) -> DVec2 {
    DVec2::new(-v.y, v.x)
}

fn orthonormal_frame(axis: DVec3, ref_dir: DVec3) -> (DVec3, DVec3, DVec3) {
    let axis = axis.normalize_or_zero();
    let mut x_axis = ref_dir - axis * ref_dir.dot(axis);
    if x_axis.length_squared() <= 1e-24 {
        x_axis = any_perpendicular(axis);
    } else {
        x_axis = x_axis.normalize();
    }
    let y_axis = axis.cross(x_axis).normalize_or_zero();
    (axis, x_axis, y_axis)
}

pub fn transform_curve(curve: &Curve3, loc: &glam::DAffine3) -> Curve3 {
    match curve {
        Curve3::Line(l) => Curve3::Line(Line3::new(loc.transform_point3(l.origin), loc.transform_vector3(l.direction))),
        Curve3::Circle(c) => {
            let center = loc.transform_point3(c.center);
            let normal = loc.transform_vector3(c.normal).normalize_or_zero();
            let x_dir = loc.transform_vector3(c.x_dir).normalize_or_zero();
            let y_dir = loc.transform_vector3(c.y_dir).normalize_or_zero();
            // OCCT-aligned: radius scales by sqrt(scale_in_plane),
            // NOT by normal length.  Use average of x_dir and y_dir
            // transform lengths (in-plane scale factors).
            let sx = loc.transform_vector3(c.x_dir).length().max(1e-12);
            let sy = loc.transform_vector3(c.y_dir).length().max(1e-12);
            let radius = c.radius * (sx * sy).sqrt(); // geometric mean = area scale
            Curve3::Circle(Circle3 {
                center,
                normal,
                x_dir,
                y_dir,
                radius,
            })
        }
        Curve3::BSpline(bs) => Curve3::BSpline(BSplineCurve3 {
            degree: bs.degree,
            knots: bs.knots.clone(),
            control_points: bs
                .control_points
                .iter()
                .map(|&p| loc.transform_point3(p))
                .collect(),
            weights: bs.weights.clone(),
            is_periodic: false,
        }),
        other => other.clone(),
    }
}

/// OCCT-aligned: Geom2d_Curve::Translate(gp_Vec2d) = Geom2d_Geometry::Translate
/// (Geom2d_Geometry.cxx L64-70, Transform with translation trsf), per-type
/// Transform implementations. Line/Circle/Ellipse/Hyperbola/Parabola translate
/// the conic position; Bezier/BSpline transform poles; TrimmedCurve translates
/// the basis and keeps the trim (Geom2d_TrimmedCurve.cxx L287-293, the
/// re-trimmed parameters are unchanged for a translation since
/// TransformedParameter(U, translation) = U); OffsetCurve translates the basis
/// and keeps the offset (Geom2d_OffsetCurve.cxx L414-419, offsetValue *= |ScaleFactor|
/// = unchanged for a translation). The rcad-only analytic types
/// (SineWave, CircleInvolute, ArchimedeanSpiral, LogarithmicSpiral, AHTBezier,
/// TBezier) have no OCCT Geom2d counterpart; they translate their defining
/// elements (SineWave is translation-invariant as it is stored phase-only).
pub fn translate_curve2d(curve: &Curve2d, offset: DVec2) -> Curve2d {
    match curve {
        Curve2d::Line(l) => Curve2d::Line(l.translate(offset)),
        Curve2d::Circle(c) => Curve2d::Circle(Circle2d {
            center: c.center + offset,
            ..*c
        }),
        Curve2d::Ellipse(e) => Curve2d::Ellipse(Ellipse2d {
            center: e.center + offset,
            ..*e
        }),
        Curve2d::Parabola(p) => Curve2d::Parabola(Parabola2d {
            origin: p.origin + offset,
            ..*p
        }),
        Curve2d::Hyperbola(h) => Curve2d::Hyperbola(Hyperbola2d {
            center: h.center + offset,
            ..*h
        }),
        Curve2d::CircleInvolute(c) => Curve2d::CircleInvolute(CircleInvolute2d {
            center: c.center + offset,
            ..*c
        }),
        Curve2d::ArchimedeanSpiral(s) => Curve2d::ArchimedeanSpiral(ArchimedeanSpiral2d {
            center: s.center + offset,
            ..*s
        }),
        Curve2d::LogarithmicSpiral(s) => Curve2d::LogarithmicSpiral(LogarithmicSpiral2d {
            center: s.center + offset,
            ..*s
        }),
        Curve2d::SineWave(s) => Curve2d::SineWave(*s),
        Curve2d::BSpline(b) => Curve2d::BSpline(BSplineCurve2 {
            control_points: b.control_points.iter().map(|&p| p + offset).collect(),
            ..b.clone()
        }),
        Curve2d::Bezier(b) => Curve2d::Bezier(BezierCurve2 {
            control_points: b.control_points.iter().map(|&p| p + offset).collect(),
            ..b.clone()
        }),
        Curve2d::Trimmed(t) => Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(translate_curve2d(&t.curve, offset)),
            ..t.clone()
        }),
        Curve2d::Offset(o) => Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(translate_curve2d(&o.basis, offset)),
            ..o.clone()
        }),
        Curve2d::AHTBezier(b) => Curve2d::AHTBezier(AHTBezierCurve2 {
            control_points: b.control_points.iter().map(|&p| p + offset).collect(),
            ..b.clone()
        }),
        Curve2d::TBezier(b) => Curve2d::TBezier(TBezierCurve2 {
            control_points: b.control_points.iter().map(|&p| p + offset).collect(),
            ..b.clone()
        }),
    }
}

/// OCCT-aligned: Geom2d_Curve::Reverse — the curve traversed in the opposite
/// direction (ReversedParameter semantics).  Used by
/// BOPTools_AlgoTools2D::AttachExistingPCurve (BOPTools_AlgoTools2D_1.cxx L80):
/// aC2DoldC->Reverse() when IsSplitToReverse.
pub fn reverse_curve2d(curve: &Curve2d) -> Curve2d {
    match curve {
        Curve2d::Line(l) => Curve2d::Line(l.with_direction(-l.direction)),
        // Geom2d_Circle::Reverse: P'(t) = P(-t) — the x frame axis is kept,
        // the y axis is negated (sin(-t) = -sin t).
        Curve2d::Circle(c) => Curve2d::Circle(Circle2d {
            y_dir: -c.y_dir,
            ..*c
        }),
        // Geom2d_Ellipse::Reverse (gp_Elips2d::Reverse): the X (major)
        // direction is kept, the stored Y direction is negated.
        Curve2d::Ellipse(e) => Curve2d::Ellipse(Ellipse2d {
            minor_dir: -e.minor_dir,
            ..*e
        }),
        Curve2d::Parabola(p) => Curve2d::Parabola(Parabola2d {
            axis_dir: -p.axis_dir,
            ..*p
        }),
        Curve2d::Hyperbola(h) => Curve2d::Hyperbola(Hyperbola2d {
            major_dir: -h.major_dir,
            ..*h
        }),
        // Geom2d_TrimmedCurve::Reverse: reverse the basis and swap the trimmed
        // bounds through ReversedParameter.
        Curve2d::Trimmed(t) => {
            let b = reverse_curve2d(&t.curve);
            let new_min = b.reversed_parameter(t.t_max);
            let new_max = b.reversed_parameter(t.t_min);
            Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(b),
                t_min: new_min,
                t_max: new_max,
            })
        }
        // Geom2d_BSplineCurve::Reverse: poles reversed, knots mirrored.
        Curve2d::BSpline(b) => {
            let mut knots: Vec<f64> = b.knots.iter().map(|k| -k).collect();
            knots.reverse();
            let k0 = knots.first().copied().unwrap_or(0.0);
            for k in knots.iter_mut() {
                *k -= k0;
            }
            Curve2d::BSpline(BSplineCurve2 {
                degree: b.degree,
                knots,
                control_points: b.control_points.iter().rev().cloned().collect(),
                weights: b.weights.iter().rev().cloned().collect(),
                // Geom2d_BSplineCurve::Reverse keeps myPeriodic.
                is_periodic: b.is_periodic,
            })
        }
        Curve2d::Bezier(b) => Curve2d::Bezier(BezierCurve2 {
            control_points: b.control_points.iter().rev().cloned().collect(),
            weights: b.weights.iter().rev().cloned().collect(),
        }),
        Curve2d::Offset(o) => Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(reverse_curve2d(&o.basis)),
            offset_distance: -o.offset_distance,
        }),
        other => other.clone(),
    }
}

/// OCCT-aligned: GeomLib::SameRange (GeomLib.cxx L842-970) for 2D curves —
/// re-parameterize `theCurve` (defined on [X1,X2]) onto [Y1,Y2] keeping the
/// geometry.  Line: translation by dU*D (L864-871); Circle: rotation (L872-889);
/// TrimmedCurve: recurse into the basis and re-trim (L890-901); other types:
/// CurveToBSplineCurve + BSplCLib::Reparametrize (L908-922, L924-969).
pub fn same_range_2d(
    tolerance: f64,
    c: Curve2d,
    x1: f64,
    x2: f64,
    y1: f64,
    y2: f64,
) -> Option<Curve2d> {
    use crate::geom::Curve2dEval;
    let tol = tolerance;
    // L854-859: ranges already equal -> the curve itself.
    if (x2 - y2).abs() <= tol && (x1 - y1).abs() <= tol {
        return Some(c);
    }
    // L862: the parameterization length must be preserved.
    let len_eq = (x2 - x1 - y2 + y1).abs() <= tol;
    if len_eq {
        match c {
            // L864-871: Line->Translate((FirstOnCurve - RequestedFirst) * D).
            Curve2d::Line(l) => {
                let d_u = x1 - y1;
                Some(Curve2d::Line(l.translate(l.direction * d_u)))
            }
            // L872-889: Circle rotation by dU (sign by IsDirect).
            Curve2d::Circle(mut cir) => {
                let is_direct =
                    (cir.x_dir.x * cir.y_dir.y - cir.x_dir.y * cir.y_dir.x) > 0.0;
                let d_u = if is_direct { x1 - y1 } else { y1 - x1 };
                cir.rotate_center(d_u);
                Some(Curve2d::Circle(cir))
            }
            // L890-901: recurse into the basis, re-trim to [Y1,Y2].
            Curve2d::Trimmed(tc) => {
                let b = same_range_2d(tolerance, *tc.curve, x1, x2, y1, y2)?;
                Some(Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(b),
                    t_min: y1,
                    t_max: y2,
                }))
            }
            // L902-922: guarded conversion — TrimmedCurve(X1,X2) -> BSpline ->
            // BSplCLib::Reparametrize(RequestedFirst, RequestedLast, Knots).
            other => {
                if (x2 - x1).abs() > crate::core::precision::PCONFUSION
                    || (y2 + y1).abs() > crate::core::precision::PCONFUSION
                {
                    let tc = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(other),
                        t_min: x1,
                        t_max: x2,
                    });
                    let mut bs = curve_to_bspline_2d(&tc)?;
                    crate::math::bspl_lib::reparametrize(y1, y2, &mut bs.knots);
                    return Some(Curve2d::BSpline(bs));
                }
                // OCCT leaves NewCurvePtr untouched (the caller's null
                // handle) when both parametric extents are below PConfusion.
                None
            }
        }
    } else {
        // L924-969: segment the curve, then BSpline + Reparametrize. The
        // PERIODIC basis keeps the requested bounds (any parameter is valid);
        // a non-periodic basis clips them to its own domain.
        let a_c_check_periodic = {
            let basis: &Curve2d = match &c {
                Curve2d::Trimmed(tc) => &tc.curve,
                other => other,
            };
            basis.is_periodic()
        };
        let [f0, f1] = c.default_domain();
        let (t_min, t_max) = if a_c_check_periodic {
            if (x2 - x1).abs() > crate::core::precision::PCONFUSION {
                (x1, x2)
            } else {
                (f0, f1)
            }
        } else {
            let u_deb = f0.max(x1);
            let u_fin = f1.min(x2);
            if (u_fin - u_deb).abs() > crate::core::precision::PCONFUSION {
                (u_deb, u_fin)
            } else {
                (f0, f1)
            }
        };
        let tc = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(c),
            t_min,
            t_max,
        });
        let mut bs = curve_to_bspline_2d(&tc)?;
        crate::math::bspl_lib::reparametrize(y1, y2, &mut bs.knots);
        Some(Curve2d::BSpline(bs))
    }
}

/// OCCT-aligned: Geom2dConvert::CurveToBSplineCurve — an exact degree-1
/// BSpline for a line, an interpolated BSpline otherwise (boolean section
/// pcurves are analytic or already BSplines; the interpolant passes through
/// the sampled points like the stored curve).
fn curve_to_bspline_2d(c: &Curve2d) -> Option<BSplineCurve2> {
    use crate::geom::Curve2dEval;
    match c {
        Curve2d::Line(l) => {
            let (t0, t1) = match c {
                Curve2d::Trimmed(tc) => (tc.t_min, tc.t_max),
                _ => (0.0, 1.0),
            };
            let p0 = l.point_at(t0);
            let p1 = l.point_at(t1);
            Some(BSplineCurve2 {
                degree: 1,
                knots: vec![t0, t0, t1, t1],
                control_points: vec![p0, p1],
                weights: vec![1.0, 1.0],
                is_periodic: false,
            })
        }
        _ => {
            let [t0, t1] = match c {
                Curve2d::Trimmed(tc) => [tc.t_min, tc.t_max],
                _ => c.default_domain(),
            };
            // OCCT Precision::IsInfinite (Precision.hxx L350-353): the
            // unbounded domains (Line2d/Parabola2d/Hyperbola2d =
            // +/-2e100) cannot be sampled as finite ranges.
            if crate::core::precision::is_infinite_value(t0)
                || crate::core::precision::is_infinite_value(t1)
                || (t1 - t0).abs() < 1e-30
            {
                return None;
            }
            let n = 23usize;
            let mut pts: Vec<glam::DVec2> = Vec::with_capacity(n + 1);
            for i in 0..=n {
                pts.push(c.point_at(t0 + (t1 - t0) * i as f64 / n as f64));
            }
            crate::base::geom_api::interpolate_points_2d(&pts).ok()
        }
    }
}

/// OCCT-aligned: apply TopLoc_Location transform to a Surface3.
pub fn transform_surface(surface: &Surface3, loc: &glam::DAffine3) -> Surface3 {
    match surface {
        Surface3::Plane(p) => Surface3::Plane(Plane::new(
            loc.transform_point3(p.origin),
            loc.transform_vector3(p.normal).normalize_or_zero(),
        )),
        Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
            origin: loc.transform_point3(c.origin),
            axis: loc.transform_vector3(c.axis).normalize_or_zero(),
            radius: c.radius * loc.transform_vector3(c.axis).length().max(1e-12),
            ref_dir: loc.transform_vector3(c.ref_dir).normalize_or_zero(),
            y_dir: c.y_dir.map(|y| loc.transform_vector3(y).normalize_or_zero()),
        }),
        Surface3::Sphere(s) => Surface3::Sphere(SphericalSurface {
            center: loc.transform_point3(s.center),
            axis: loc.transform_vector3(s.axis).normalize_or_zero(),
            radius: s.radius * loc.transform_vector3(s.axis).length().max(1e-12),
            ref_dir: loc.transform_vector3(s.ref_dir).normalize_or_zero(),
        }),
        Surface3::Cone(c) => Surface3::Cone(ConicalSurface {
            apex: loc.transform_point3(c.apex),
            axis: loc.transform_vector3(c.axis).normalize_or_zero(),
            radius: c.radius * loc.transform_vector3(c.axis).length().max(1e-12),
            half_angle_rad: c.half_angle_rad,
            ref_dir: loc.transform_vector3(c.ref_dir).normalize_or_zero(),
        }),
        Surface3::Torus(t) => Surface3::Torus(ToroidalSurface {
            center: loc.transform_point3(t.center),
            axis: loc.transform_vector3(t.axis).normalize_or_zero(),
            ref_dir: loc.transform_vector3(t.ref_dir).normalize_or_zero(),
            major_radius: t.major_radius * loc.transform_vector3(t.axis).length().max(1e-12),
            minor_radius: t.minor_radius,
        }),
        Surface3::BSpline(bs) => Surface3::BSpline(BSplineSurface {
            degree_u: bs.degree_u,
            degree_v: bs.degree_v,
            knots_u: bs.knots_u.clone(),
            knots_v: bs.knots_v.clone(),
            control_points: bs
                .control_points
                .iter()
                .map(|row| row.iter().map(|&p| loc.transform_point3(p)).collect())
                .collect(),
            weights: bs.weights.clone(),
            is_periodic_u: bs.is_periodic_u,
            is_periodic_v: bs.is_periodic_v,
        }),
        other => other.clone(),
    }
}

pub mod bspline2d_dn;
pub mod bspline_ops;
pub mod bspline_surface_ops;
pub mod curve_dn;
pub mod eval;
pub mod eval_b;
pub mod eval_c;
pub mod extrusion_utils;
pub mod offset2d_dn;
pub mod offset_surface_utils;
pub mod offset_surface_utils_b;
pub mod osculating_surface;
pub mod revolution_utils;
#[cfg(test)]
pub mod tests;
