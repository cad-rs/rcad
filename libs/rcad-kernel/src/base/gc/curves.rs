//! 3D curve construction algorithms.
//!
//! OCCT GC package: GC_MakeCircle, GC_MakeLine, GC_MakeSegment, GC_MakeEllipse,
//! GC_MakeHyperbola, GC_MakeArcOfCircle, GC_MakeArcOfEllipse, GC_MakeArcOfHyperbola,
//! GC_MakeArcOfParabola.

#![allow(clippy::manual_clamp)]

use crate::geom::{
    Circle3, Curve3, Ellipse3, Hyperbola3, Line3, Parabola3, Point3, TrimmedCurve3, Vec3,
};

use super::{GceError, TOL_CONF, TOL_FLOAT, points_coincident};

// ============================================================================
// GC_MakeCircle
// ============================================================================

/// Construct a circle from center, normal direction, and radius.
///
/// OCCT: `GC_MakeCircle(gp_Pnt, gp_Dir, double)`.
/// Returns `GceError::NegativeRadius` when radius < 0.
pub fn make_circle(center: Point3, normal: Vec3, radius: f64) -> Result<Circle3, GceError> {
    if radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    Ok(Circle3::new(center, normal, radius))
}

/// Construct a circle from three points.
///
/// OCCT: `GC_MakeCircle(gp_Pnt, gp_Pnt, gp_Pnt)`.
/// The three points define the circle (they must not be collinear or coincident).
/// Returns `GceError::ColinearPoints` or `GceError::ConfusedPoints` on invalid input.
pub fn make_circle_3p(p1: Point3, p2: Point3, p3: Point3) -> Result<Circle3, GceError> {
    // OCCT: gce_MakeCirc(P1, P2, P3)
    if points_coincident(p1, p2) || points_coincident(p1, p3) || points_coincident(p2, p3) {
        return Err(GceError::ConfusedPoints);
    }

    // Midpoints of chords
    let m1 = (p1 + p2) * 0.5;
    let m2 = (p2 + p3) * 0.5;

    // Direction vectors of chords
    let d1 = p2 - p1;
    let d2 = p3 - p2;

    // Normal to the plane of the three points
    let normal = d1.cross(d2);
    if normal.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ColinearPoints);
    }
    let normal = normal.normalize();

    // Perpendicular bisectors of chords
    let perp1 = normal.cross(d1).normalize();
    let perp2 = normal.cross(d2).normalize();

    // Find center by intersecting the two perpendicular bisectors
    // Solve: m1 + t1 * perp1 = m2 + t2 * perp2
    // => t1 * perp1 - t2 * perp2 = m2 - m1
    let diff = m2 - m1;
    let cross_p1_p2 = perp1.cross(perp2);
    if cross_p1_p2.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ColinearPoints);
    }
    let t1 = diff.cross(perp2).dot(cross_p1_p2) / cross_p1_p2.length_squared();
    let center = m1 + perp1 * t1;

    let radius = (p1 - center).length();
    if radius < TOL_CONF {
        return Err(GceError::NullRadius);
    }

    Ok(Circle3::new(center, normal, radius))
}

/// Construct a circle from center, axis point (defines normal), and radius.
///
/// OCCT: `GC_MakeCircle(gp_Pnt, gp_Pnt, double)`.
/// The normal direction is defined by vector from center to axis_point.
pub fn make_circle_center_axis(center: Point3, axis_point: Point3, radius: f64) -> Result<Circle3, GceError> {
    if radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let normal = axis_point - center;
    if normal.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    Ok(Circle3::new(center, normal, radius))
}

/// Construct a circle from an axis line (location + direction) and radius.
///
/// OCCT: `GC_MakeCircle(gp_Ax1, double)`.
/// The axis location is the center; the axis direction is the normal.
pub fn make_circle_axis(center: Point3, axis_dir: Vec3, radius: f64) -> Result<Circle3, GceError> {
    if radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let normal = axis_dir.normalize_or_zero();
    if normal.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Circle3::new(center, normal, radius))
}

/// Construct a circle concentric to another circle with a radius offset.
///
/// OCCT: `GC_MakeCircle(gp_Circ, double)`.
pub fn make_circle_concentric_offset(circle: &Circle3, dist: f64) -> Result<Circle3, GceError> {
    let new_radius = circle.radius + dist;
    if new_radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    Ok(Circle3 {
        center: circle.center,
        normal: circle.normal,
        x_dir: circle.x_dir,
        y_dir: circle.y_dir,
        radius: new_radius,
    })
}

/// Construct a circle concentric to another circle passing through a point.
///
/// OCCT: `GC_MakeCircle(gp_Circ, gp_Pnt)`.
pub fn make_circle_concentric_point(circle: &Circle3, point: Point3) -> Result<Circle3, GceError> {
    let d = point - circle.center;
    let planar = d - circle.normal * d.dot(circle.normal);
    let radius = planar.length();
    if radius < TOL_CONF {
        return Err(GceError::NullRadius);
    }
    Ok(Circle3 {
        center: circle.center,
        normal: circle.normal,
        x_dir: circle.x_dir,
        y_dir: circle.y_dir,
        radius,
    })
}

// ============================================================================
// GC_MakeLine
// ============================================================================

/// Construct a line from an axis (location + direction).
///
/// OCCT: `GC_MakeLine(gp_Ax1)`.
pub fn make_line_axis(origin: Point3, direction: Vec3) -> Result<Line3, GceError> {
    let dir = direction.normalize_or_zero();
    if dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Line3 { origin, direction: dir })
}

/// Construct a line from point and direction.
///
/// OCCT: `GC_MakeLine(gp_Pnt, gp_Dir)`.
pub fn make_line_pd(point: Point3, direction: Vec3) -> Result<Line3, GceError> {
    let dir = direction.normalize_or_zero();
    if dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Line3 { origin: point, direction: dir })
}

/// Construct a line passing through two points.
///
/// OCCT: `GC_MakeLine(gp_Pnt, gp_Pnt)`.
/// Returns `GceError::ConfusedPoints` when the two points are coincident.
pub fn make_line_2p(p1: Point3, p2: Point3) -> Result<Line3, GceError> {
    let direction = p2 - p1;
    if direction.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    Ok(Line3::new(p1, direction))
}

/// Construct a line parallel to an existing line passing through a point.
///
/// OCCT: `GC_MakeLine(gp_Lin, gp_Pnt)`.
pub fn make_line_parallel_point(line: &Line3, point: Point3) -> Result<Line3, GceError> {
    Ok(Line3 {
        origin: point,
        direction: line.direction,
    })
}

// ============================================================================
// GC_MakeSegment — constructs a trimmed line (Curve3::Trimmed(Line3))
// ============================================================================

/// Construct a segment from two points (trimmed line between them).
///
/// OCCT: `GC_MakeSegment(gp_Pnt, gp_Pnt)`.
pub fn make_segment(p1: Point3, p2: Point3) -> Result<Curve3, GceError> {
    let direction = p2 - p1;
    if direction.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    let line = Line3::new(p1, direction);
    // Parameter values at the endpoints
    let t1 = 0.0;
    let t2 = (p2 - p1).length();
    Ok(Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Line(line),
        t1,
        t2,
    )))
}

/// Construct a segment from a line and two parameter values.
///
/// OCCT: `GC_MakeSegment(gp_Lin, double, double)`.
pub fn make_segment_on_line(line: &Line3, p1: f64, p2: f64) -> Result<Curve3, GceError> {
    if (p1 - p2).abs() < TOL_FLOAT {
        return Err(GceError::NullLength);
    }
    let (t1, t2) = if p1 < p2 { (p1, p2) } else { (p2, p1) };
    Ok(Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Line(*line),
        t1,
        t2,
    )))
}

// ============================================================================
// GC_MakeEllipse
// ============================================================================

/// Construct an ellipse from center, normal, major direction, and radii.
///
/// OCCT: `GC_MakeEllipse(gp_Ax2, double, double)`.
pub fn make_ellipse(center: Point3, normal: Vec3, major_dir: Vec3, major_radius: f64, minor_radius: f64) -> Result<Ellipse3, GceError> {
    if major_radius < 0.0 || minor_radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let normal = normal.normalize_or_zero();
    if normal.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    let major_dir = major_dir.normalize_or_zero();
    if major_dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Ellipse3 {
        center,
        normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}

// ============================================================================
// GC_MakeHyperbola
// ============================================================================

/// Construct a hyperbola from center, normal, major direction, and semi-axes.
///
/// OCCT: `GC_MakeHyperbola(gp_Ax2, double, double)`.
pub fn make_hyperbola(center: Point3, normal: Vec3, major_dir: Vec3, semi_major: f64, semi_minor: f64) -> Result<Hyperbola3, GceError> {
    if semi_major < 0.0 || semi_minor < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let normal = normal.normalize_or_zero();
    if normal.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    let major_dir = major_dir.normalize_or_zero();
    if major_dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Hyperbola3 {
        center,
        normal,
        major_dir,
        semi_major,
        semi_minor,
    })
}

// ============================================================================
// GC_MakeArcOfCircle — constructs a Curve3::Trimmed containing a Circle3
// ============================================================================

/// Construct a circular arc from three points (start, intermediate, end).
///
/// OCCT: `GC_MakeArcOfCircle(gp_Pnt, gp_Pnt, gp_Pnt)`.
pub fn make_arc_of_circle_3p(p1: Point3, p2: Point3, p3: Point3) -> Result<Curve3, GceError> {
    // Find the circle through the three points
    let circle = make_circle_3p(p1, p2, p3)?;

    // Compute parameters on the circle
    let center = circle.center;
    let x_dir = circle.x_dir;
    let y_dir = circle.y_dir;

    fn param_on_circle(p: Point3, center: Point3, x_dir: Vec3, y_dir: Vec3) -> f64 {
        let d = p - center;
        let x = d.dot(x_dir);
        let y = d.dot(y_dir);
        y.atan2(x)
    }

    let t_start = param_on_circle(p1, center, x_dir, y_dir);
    let t_end = param_on_circle(p3, center, x_dir, y_dir);
    let t_mid = param_on_circle(p2, center, x_dir, y_dir);

    // OCCT: ensure the arc goes from p1 to p3 passing through p2 (mid parameter order)
    // If p2 is not between p1 and p3 in the natural direction, adjust period.
    let two_pi = 2.0 * std::f64::consts::PI;
    let (mut start, mut end) = (t_start, t_end);

    // Adjust so mid lies between start and end
    let in_range = |t: f64, a: f64, b: f64| -> bool {
        if a <= b {
            t >= a && t <= b
        } else {
            t >= a || t <= b
        }
    };

    if !in_range(t_mid, t_start, t_end) {
        if t_end > t_start {
            end = t_end - two_pi;
        } else {
            start = t_start - two_pi;
            end = t_end;
        }
    }

    if (start - end).abs() < TOL_FLOAT {
        return Err(GceError::ConfusedPoints);
    }

    Ok(Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Circle(circle),
        start.min(end),
        start.max(end),
    )))
}

// Note: t3 is not defined in the above code, fix in next edit
// Actually the variable should be t3 vs t2 confusion. Let me fix.
// OCCT: p1=start, p2=intermediate, p3=end.
// So t1 = param(p1), t3 = param(p3), t_mid = param(p2).

/// Construct a circular arc from a circle and two parameter values.
///
/// OCCT: `GC_MakeArcOfCircle(gp_Circ, double, double, bool)`.
pub fn make_arc_of_circle_params(circle: &Circle3, alpha1: f64, alpha2: f64) -> Result<Curve3, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    let (t1, t2) = if alpha1 < alpha2 { (alpha1, alpha2) } else { (alpha2, alpha1) };
    Ok(Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Circle(*circle),
        t1,
        t2,
    )))
}

/// Construct a circular arc from two points and a center.
///
/// The center defines the circle; p1 and p2 define the arc endpoints.
/// OCCT: `GC_MakeArcOfCircle(gp_Pnt, gp_Pnt, gp_Pnt)` — center variant.
/// Note: this follows OCCT's gp_Circ from center + 2 points overload.
pub fn make_arc_of_circle_center(start: Point3, end: Point3, center: Point3) -> Result<Curve3, GceError> {
    let d1 = start - center;
    let d2 = end - center;
    let r1 = d1.length();
    let r2 = d2.length();
    if r1 < TOL_CONF || r2 < TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    if (r1 - r2).abs() > TOL_CONF {
        return Err(GceError::ConfusedPoints); // Not on same circle
    }

    let radius = (r1 + r2) * 0.5;
    let normal = d1.cross(d2);
    if normal.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ColinearPoints);
    }

    let circle = Circle3::new(center, normal, radius);

    let x_dir = circle.x_dir;
    let y_dir = circle.y_dir;
    let t1 = (d1.dot(y_dir)).atan2(d1.dot(x_dir));
    let t2 = (d2.dot(y_dir)).atan2(d2.dot(x_dir));

    Ok(Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Circle(circle),
        t1.min(t2),
        t1.max(t2),
    )))
}

// ============================================================================
// GC_MakeArcOfEllipse
// ============================================================================

/// Construct an elliptical arc on an ellipse between two parameter values.
///
/// OCCT: `GC_MakeArcOfEllipse(gp_Elips, double, double, bool)`.
pub fn make_arc_of_ellipse(ellipse: &Ellipse3, alpha1: f64, alpha2: f64) -> Result<Curve3, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    // Wrap ellipse in a trimmed curve
    let ellipse_curve = Curve3::Ellipse(*ellipse);
    // The parameter range for an ellipse in CurveEval is [0, 2*pi]; clamp
    let two_pi = 2.0 * std::f64::consts::PI;
    let t1 = alpha1.rem_euclid(two_pi);
    let t2 = alpha2.rem_euclid(two_pi);
    let (t_min, t_max) = if t1 <= t2 { (t1, t2) } else { (t2, t1) };
    Ok(Curve3::Trimmed(TrimmedCurve3::new(ellipse_curve, t_min, t_max)))
}

// ============================================================================
// GC_MakeArcOfHyperbola
// ============================================================================

/// Construct an arc of hyperbola between two parameter values.
///
/// OCCT: `GC_MakeArcOfHyperbola(gp_Hypr, double, double, bool)`.
pub fn make_arc_of_hyperbola(hyperbola: &Hyperbola3, alpha1: f64, alpha2: f64) -> Result<Curve3, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    let (t1, t2) = if alpha1 < alpha2 { (alpha1, alpha2) } else { (alpha2, alpha1) };
    Ok(Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Hyperbola(*hyperbola),
        t1,
        t2,
    )))
}

// ============================================================================
// GC_MakeArcOfParabola
// ============================================================================

/// Construct an arc of parabola between two parameter values.
///
/// OCCT: `GC_MakeArcOfParabola(gp_Parab, double, double, bool)`.
/// Note: OCCT only has `GC_MakeArcOfParabola2d` — the 3D variant wraps
/// `Geom_Parabola` with a parameter range.
pub fn make_arc_of_parabola(parabola: &Parabola3, alpha1: f64, alpha2: f64) -> Result<Curve3, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    let (t1, t2) = if alpha1 < alpha2 { (alpha1, alpha2) } else { (alpha2, alpha1) };
    Ok(Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Parabola(*parabola),
        t1,
        t2,
    )))
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::CurveEval;

    #[test]
    fn test_make_circle_basic() {
        let c = make_circle(DVec3::ZERO, DVec3::Z, 5.0).unwrap();
        assert!((c.radius - 5.0).abs() < 1e-12);
        assert!((c.normal - DVec3::Z).length() < 1e-12);
    }

    #[test]
    fn test_make_circle_negative_radius() {
        assert_eq!(
            make_circle(DVec3::ZERO, DVec3::Z, -1.0).unwrap_err(),
            GceError::NegativeRadius
        );
    }

    #[test]
    fn test_make_circle_3p() {
        let c = make_circle_3p(
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(-1.0, 0.0, 0.0),
        )
        .unwrap();
        assert!((c.radius - 1.0).abs() < 1e-7);
        assert!((c.center - DVec3::ZERO).length() < 1e-7);
    }

    #[test]
    fn test_make_circle_3p_collinear() {
        assert_eq!(
            make_circle_3p(
                DVec3::ZERO,
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(2.0, 0.0, 0.0),
            )
            .unwrap_err(),
            GceError::ColinearPoints
        );
    }

    #[test]
    fn test_make_line_2p() {
        let l = make_line_2p(DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)).unwrap();
        assert!((l.direction - DVec3::X).length() < 1e-12);
    }

    #[test]
    fn test_make_line_2p_confused() {
        assert_eq!(
            make_line_2p(DVec3::new(1.0, 2.0, 3.0), DVec3::new(1.0, 2.0, 3.0)).unwrap_err(),
            GceError::ConfusedPoints
        );
    }

    #[test]
    fn test_make_line_pd() {
        let l = make_line_pd(DVec3::new(1.0, 2.0, 3.0), DVec3::new(0.0, 1.0, 0.0)).unwrap();
        assert!((l.direction - DVec3::Y).length() < 1e-12);
        assert!((l.origin - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_line_null_axis() {
        assert_eq!(
            make_line_pd(DVec3::ZERO, DVec3::ZERO).unwrap_err(),
            GceError::NullAxis
        );
    }

    #[test]
    fn test_make_segment() {
        let seg = make_segment(DVec3::ZERO, DVec3::new(3.0, 0.0, 0.0)).unwrap();
        if let Curve3::Trimmed(tc) = &seg {
            assert!((tc.first - 0.0).abs() < 1e-12);
            assert!((tc.last - 3.0).abs() < 1e-12);
        } else {
            panic!("expected Trimmed curve");
        }
        // Endpoints evaluate correctly
        let p0 = seg.point_at(0.0);
        let p3 = seg.point_at(3.0);
        assert!((p0 - DVec3::ZERO).length() < 1e-12);
        assert!((p3 - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_ellipse() {
        let e = make_ellipse(
            DVec3::ZERO,
            DVec3::Z,
            DVec3::X,
            5.0,
            3.0,
        )
        .unwrap();
        assert!((e.major_radius - 5.0).abs() < 1e-12);
        assert!((e.minor_radius - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_make_hyperbola() {
        let h = make_hyperbola(
            DVec3::ZERO,
            DVec3::Z,
            DVec3::X,
            4.0,
            2.0,
        )
        .unwrap();
        assert!((h.semi_major - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_make_arc_of_circle_params() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let arc = make_arc_of_circle_params(&circle, 0.0, std::f64::consts::PI).unwrap();
        if let Curve3::Trimmed(tc) = &arc {
            assert!((tc.first - 0.0).abs() < 1e-12);
            assert!((tc.last - std::f64::consts::PI).abs() < 1e-12);
        } else {
            panic!("expected Trimmed curve");
        }
    }

    #[test]
    fn test_make_arc_of_ellipse() {
        let ellipse = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 4.0,
            minor_radius: 2.0,
        };
        let arc = make_arc_of_ellipse(&ellipse, 0.0, std::f64::consts::PI * 0.5).unwrap();
        if let Curve3::Trimmed(tc) = &arc {
            assert!((tc.first - 0.0).abs() < 1e-12);
        } else {
            panic!("expected Trimmed curve");
        }
    }

    #[test]
    fn test_make_arc_of_hyperbola() {
        let h = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };
        let arc = make_arc_of_hyperbola(&h, 0.0, 1.0).unwrap();
        if let Curve3::Trimmed(_) = &arc {
            // ok
        } else {
            panic!("expected Trimmed curve");
        }
    }

    #[test]
    fn test_make_circle_concentric_point() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 3.0);
        let pt = DVec3::new(0.0, 5.0, 0.0);
        let c = make_circle_concentric_point(&circle, pt).unwrap();
        assert!((c.radius - 5.0).abs() < 1e-12);
        assert!((c.center - circle.center).length() < 1e-12);
    }
}
