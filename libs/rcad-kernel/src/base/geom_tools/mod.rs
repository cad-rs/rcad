//! Curve and surface dump/read utilities (GeomTools).
//!
//! OCCT TKGeomBase GeomTools package.
//! Provides text dump and debug-print for geometry types.

use crate::geom::{Curve3, Surface3, Curve2d};

/// Dump a 3D curve to a string (human-readable).
///
/// OCCT: `GeomTools::Dump(Curve)`.
pub fn dump_curve(curve: &Curve3) -> String {
    match curve {
        Curve3::Line(l) => format!("Line(origin: {:?}, dir: {:?})", l.origin, l.direction),
        Curve3::Circle(c) => format!("Circle(center: {:?}, normal: {:?}, r: {})", c.center, c.normal, c.radius),
        Curve3::Ellipse(e) => format!("Ellipse(center: {:?}, major_radius: {}, minor_radius: {})", e.center, e.major_radius, e.minor_radius),
        Curve3::BSpline(b) => format!("BSpline(degree: {}, poles: {})", b.degree, b.control_points.len()),
        Curve3::Bezier(b) => format!("Bezier(poles: {})", b.control_points.len()),
        Curve3::Trimmed(tc) => format!("Trimmed(basis: {}, range: [{}, {}])", dump_curve(&tc.curve), tc.first, tc.last),
        Curve3::Parabola(p) => format!("Parabola(vertex: {:?}, focal: {})", p.vertex, p.focal_param),
        Curve3::Hyperbola(h) => format!("Hyperbola(center: {:?}, a: {}, b: {})", h.center, h.semi_major, h.semi_minor),
        Curve3::CircularHelix(h) => format!("Helix(radius: {}, pitch: {})", h.radius, h.pitch),
        Curve3::SineWave(s) => format!("SineWave(amp: {}, freq: {})", s.amplitude, s.frequency),
        Curve3::Offset(o) => format!("OffsetCurve(dist: {})", o.offset_distance),
    }
}

/// Dump a surface to a string.
///
/// OCCT: `GeomTools::Dump(Surface)`.
pub fn dump_surface(surface: &Surface3) -> String {
    match surface {
        Surface3::Plane(p) => format!("Plane(origin: {:?}, normal: {:?})", p.origin, p.normal),
        Surface3::Cylinder(c) => format!("Cylinder(origin: {:?}, axis: {:?}, r: {})", c.origin, c.axis, c.radius),
        Surface3::Sphere(s) => format!("Sphere(center: {:?}, r: {})", s.center, s.radius),
        Surface3::Cone(c) => format!("Cone(apex: {:?}, half_angle: {})", c.apex, c.half_angle_rad),
        Surface3::Torus(t) => format!("Torus(center: {:?}, major_r: {}, minor_r: {})", t.center, t.major_radius, t.minor_radius),
        Surface3::BSpline(b) => format!("BSplineSurface(deg_u: {}, deg_v: {}, poles: {}x{})", b.degree_u, b.degree_v, b.control_points.len(), b.control_points.first().map(|r| r.len()).unwrap_or(0)),
        Surface3::Bezier(b) => format!("BezierSurface(poles: {}x{})", b.control_points.len(), b.control_points.first().map(|r| r.len()).unwrap_or(0)),
        Surface3::Trimmed(t) => format!("TrimmedSurface(trim: {:?})", t.trim),
        _ => format!("{:?}", surface),
    }
}

/// Dump a 2D curve to a string.
///
/// OCCT: `GeomTools::Dump(Curve2d)`.
pub fn dump_curve2d(curve: &Curve2d) -> String {
    match curve {
        Curve2d::Line(l) => format!("Line2d(origin: {:?}, dir: {:?})", l.origin, l.direction),
        Curve2d::Circle(c) => format!("Circle2d(center: {:?}, r: {})", c.center, c.radius),
        Curve2d::Trimmed(t) => format!("Trimmed2d(range: [{}, {}])", t.t_min, t.t_max),
        Curve2d::BSpline(b) => format!("BSpline2d(degree: {}, poles: {})", b.degree, b.control_points.len()),
        Curve2d::Ellipse(e) => format!("Ellipse2d(center: {:?}, major_r: {}, minor_r: {})", e.center, e.major_radius, e.minor_radius),
        _ => format!("{:?}", curve),
    }
}
