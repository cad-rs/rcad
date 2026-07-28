//! 2D curve construction algorithms.
//!
//! OCCT GC package: GC_MakeCircle2d, GC_MakeLine2d, GC_MakeSegment2d,
//! GC_MakeEllipse2d, GC_MakeHyperbola2d, GC_MakeParabola2d,
//! GC_MakeArcOfCircle2d, GC_MakeArcOfEllipse2d, GC_MakeArcOfHyperbola2d,
//! GC_MakeArcOfParabola2d.

#![allow(clippy::manual_clamp)]

use glam::DVec2;

use crate::geom::{
    Circle2d, Curve2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d, Point2, TrimmedCurve2,
    Vec2,
};

use super::{
    GceError, TOL_CONF, TOL_FLOAT, points_coincident_2d, vectors_parallel_2d,
};

// ============================================================================
// GC_MakeCircle2d
// ============================================================================

/// Construct a 2D circle from center and radius.
///
/// OCCT: `GC_MakeCircle2d(gp_Circ2d)`.
pub fn make_circle2d(center: Point2, radius: f64) -> Result<Circle2d, GceError> {
    if radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    Ok(Circle2d::new(center, radius))
}

/// Construct a 2D circle from three points.
///
/// OCCT: `GC_MakeCircle2d(gp_Pnt2d, gp_Pnt2d, gp_Pnt2d)`.
pub fn make_circle2d_3p(p1: Point2, p2: Point2, p3: Point2) -> Result<Circle2d, GceError> {
    if points_coincident_2d(p1, p2)
        || points_coincident_2d(p1, p3)
        || points_coincident_2d(p2, p3)
    {
        return Err(GceError::ConfusedPoints);
    }

    // Perpendicular bisector intersection in 2D
    let m1 = (p1 + p2) * 0.5;
    let m2 = (p2 + p3) * 0.5;
    let d1 = p2 - p1;
    let d2 = p3 - p2;

    // 2D perpendicular vectors (rotate 90° CCW)
    let perp1 = DVec2::new(-d1.y, d1.x);
    let perp2 = DVec2::new(-d2.y, d2.x);

    // Check collinearity
    if vectors_parallel_2d(d1, d2) {
        return Err(GceError::ColinearPoints);
    }

    // Solve: m1 + t1 * perp1 = m2 + t2 * perp2
    // => t1 * perp1.x - t2 * perp2.x = m2.x - m1.x
    // => t1 * perp1.y - t2 * perp2.y = m2.y - m1.y
    let diff = m2 - m1;
    let det = perp1.x * perp2.y - perp1.y * perp2.x;
    if det.abs() < TOL_CONF * TOL_CONF {
        return Err(GceError::ColinearPoints);
    }
    let t1 = (diff.x * perp2.y - diff.y * perp2.x) / det;
    let center = m1 + perp1 * t1;

    let radius = (p1 - center).length();
    if radius < TOL_CONF {
        return Err(GceError::NullRadius);
    }

    Ok(Circle2d::new(center, radius))
}

/// Construct a 2D circle concentric to another circle with radius offset.
///
/// OCCT: `GC_MakeCircle2d(gp_Circ2d, double)`.
pub fn make_circle2d_concentric_offset(circle: &Circle2d, dist: f64) -> Result<Circle2d, GceError> {
    let new_radius = circle.radius + dist;
    if new_radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    Ok(Circle2d {
        center: circle.center,
        x_dir: circle.x_dir,
        y_dir: circle.y_dir,
        radius: new_radius,
    })
}

/// Construct a 2D circle concentric to another passing through a point.
///
/// OCCT: `GC_MakeCircle2d(gp_Circ2d, gp_Pnt2d)`.
pub fn make_circle2d_concentric_point(circle: &Circle2d, point: Point2) -> Result<Circle2d, GceError> {
    let radius = (point - circle.center).length();
    if radius < TOL_CONF {
        return Err(GceError::NullRadius);
    }
    Ok(Circle2d {
        center: circle.center,
        x_dir: circle.x_dir,
        y_dir: circle.y_dir,
        radius,
    })
}

// ============================================================================
// GC_MakeLine2d
// ============================================================================

/// Construct a 2D line from point and direction.
///
/// OCCT: `GC_MakeLine2d(gp_Pnt2d, gp_Dir2d)`.
pub fn make_line2d_pd(point: Point2, direction: Vec2) -> Result<Line2d, GceError> {
    let dir = direction.normalize_or_zero();
    if dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Line2d::new(point, dir))
}

/// Construct a 2D line passing through two points.
///
/// OCCT: `GC_MakeLine2d(gp_Pnt2d, gp_Pnt2d)`.
pub fn make_line2d_2p(p1: Point2, p2: Point2) -> Result<Line2d, GceError> {
    let direction = p2 - p1;
    if direction.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    Ok(Line2d::new(p1, direction))
}

/// Construct a 2D line parallel to another passing through a point.
///
/// OCCT: `GC_MakeLine2d(gp_Lin2d, gp_Pnt2d)`.
pub fn make_line2d_parallel_point(line: &Line2d, point: Point2) -> Result<Line2d, GceError> {
    Ok(Line2d {
        origin: point,
        direction: line.direction,
    })
}

// ============================================================================
// GC_MakeSegment2d
// ============================================================================

/// Construct a 2D segment from two points.
///
/// OCCT: `GC_MakeSegment2d(gp_Pnt2d, gp_Pnt2d)`.
pub fn make_segment2d(p1: Point2, p2: Point2) -> Result<Curve2d, GceError> {
    let direction = p2 - p1;
    if direction.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    let line = Line2d::new(p1, direction);
    let t1 = 0.0;
    let t2 = direction.length();
    Ok(Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(line)),
        t_min: t1,
        t_max: t2,
    }))
}

/// Construct a 2D segment from a line and two parameter values.
///
/// OCCT: `GC_MakeSegment2d(gp_Lin2d, double, double)`.
pub fn make_segment2d_on_line(line: &Line2d, p1: f64, p2: f64) -> Result<Curve2d, GceError> {
    if (p1 - p2).abs() < TOL_FLOAT {
        return Err(GceError::NullLength);
    }
    let (t1, t2) = if p1 < p2 { (p1, p2) } else { (p2, p1) };
    Ok(Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(*line)),
        t_min: t1,
        t_max: t2,
    }))
}

// ============================================================================
// GC_MakeEllipse2d
// ============================================================================

/// Construct a 2D ellipse from center, major direction, and radii.
///
/// OCCT: `GC_MakeEllipse2d(gp_Ax2d, double, double)`.
pub fn make_ellipse2d(
    center: Point2,
    major_dir: Vec2,
    major_radius: f64,
    minor_radius: f64,
) -> Result<Ellipse2d, GceError> {
    if major_radius < 0.0 || minor_radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let major_dir = major_dir.normalize_or_zero();
    if major_dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Ellipse2d {
        center,
        major_dir,
        major_radius,
        minor_radius,
    })
}

// ============================================================================
// GC_MakeHyperbola2d
// ============================================================================

/// Construct a 2D hyperbola from center, major direction, and semi-axes.
///
/// OCCT: `GC_MakeHyperbola2d(gp_Ax2d, double, double)`.
pub fn make_hyperbola2d(
    center: Point2,
    major_dir: Vec2,
    semi_major: f64,
    semi_minor: f64,
) -> Result<Hyperbola2d, GceError> {
    if semi_major < 0.0 || semi_minor < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let major_dir = major_dir.normalize_or_zero();
    if major_dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Hyperbola2d {
        center,
        major_dir,
        semi_major,
        semi_minor,
    })
}

// ============================================================================
// GC_MakeParabola2d
// ============================================================================

/// Construct a 2D parabola from origin, axis direction, and focal parameter.
///
/// OCCT: `GC_MakeParabola2d(gp_Ax2d, double)`.
pub fn make_parabola2d(origin: Point2, axis_dir: Vec2, focal_param: f64) -> Result<Parabola2d, GceError> {
    if focal_param < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let axis_dir = axis_dir.normalize_or_zero();
    if axis_dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Parabola2d {
        origin,
        axis_dir,
        focal_param,
    })
}

// ============================================================================
// GC_MakeArcOfCircle2d
// ============================================================================

/// Construct a 2D circular arc from three points (start, intermediate, end).
///
/// OCCT: `GC_MakeArcOfCircle2d(gp_Pnt2d, gp_Pnt2d, gp_Pnt2d)`.
pub fn make_arc_of_circle2d_3p(p1: Point2, p2: Point2, p3: Point2) -> Result<Curve2d, GceError> {
    let circle = make_circle2d_3p(p1, p2, p3)?;

    fn param_on_circle(p: Point2, center: Point2, x_dir: Vec2, y_dir: Vec2) -> f64 {
        let d = p - center;
        d.dot(y_dir).atan2(d.dot(x_dir))
    }

    let t_start = param_on_circle(p1, circle.center, circle.x_dir, circle.y_dir);
    let t_end = param_on_circle(p3, circle.center, circle.x_dir, circle.y_dir);
    let t_mid = param_on_circle(p2, circle.center, circle.x_dir, circle.y_dir);

    let two_pi = 2.0 * std::f64::consts::PI;
    let (mut start, mut end) = (t_start, t_end);

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

    let (t_min, t_max) = if start < end { (start, end) } else { (end, start) };
    Ok(Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Circle(circle)),
        t_min,
        t_max,
    }))
}

/// Construct a 2D circular arc from a circle and two parameter values.
///
/// OCCT: `GC_MakeArcOfCircle2d(gp_Circ2d, double, double, bool)`.
pub fn make_arc_of_circle2d_params(circle: &Circle2d, alpha1: f64, alpha2: f64) -> Result<Curve2d, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    let (t1, t2) = if alpha1 < alpha2 { (alpha1, alpha2) } else { (alpha2, alpha1) };
    Ok(Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Circle(*circle)),
        t_min: t1,
        t_max: t2,
    }))
}

// ============================================================================
// GC_MakeArcOfEllipse2d
// ============================================================================

/// Construct a 2D elliptical arc between two parameter values.
///
/// OCCT: `GC_MakeArcOfEllipse2d(gp_Elips2d, double, double, bool)`.
pub fn make_arc_of_ellipse2d(ellipse: &Ellipse2d, alpha1: f64, alpha2: f64) -> Result<Curve2d, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    let two_pi = 2.0 * std::f64::consts::PI;
    let t1 = alpha1.rem_euclid(two_pi);
    let t2 = alpha2.rem_euclid(two_pi);
    let (t_min, t_max) = if t1 <= t2 { (t1, t2) } else { (t2, t1) };
    Ok(Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Ellipse(*ellipse)),
        t_min,
        t_max,
    }))
}

// ============================================================================
// GC_MakeArcOfHyperbola2d
// ============================================================================

/// Construct a 2D arc of hyperbola between two parameter values.
///
/// OCCT: `GC_MakeArcOfHyperbola2d(gp_Hypr2d, double, double, bool)`.
pub fn make_arc_of_hyperbola2d(hyperbola: &Hyperbola2d, alpha1: f64, alpha2: f64) -> Result<Curve2d, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    let (t1, t2) = if alpha1 < alpha2 { (alpha1, alpha2) } else { (alpha2, alpha1) };
    Ok(Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Hyperbola(*hyperbola)),
        t_min: t1,
        t_max: t2,
    }))
}

// ============================================================================
// GC_MakeArcOfParabola2d
// ============================================================================

/// Construct a 2D arc of parabola between two parameter values.
///
/// OCCT: `GC_MakeArcOfParabola2d(gp_Parab2d, double, double, bool)`.
pub fn make_arc_of_parabola2d(parabola: &Parabola2d, alpha1: f64, alpha2: f64) -> Result<Curve2d, GceError> {
    if (alpha1 - alpha2).abs() < TOL_FLOAT {
        return Err(GceError::SameParameters);
    }
    let (t1, t2) = if alpha1 < alpha2 { (alpha1, alpha2) } else { (alpha2, alpha1) };
    Ok(Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Parabola(*parabola)),
        t_min: t1,
        t_max: t2,
    }))
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Curve2dEval;

    #[test]
    fn test_make_circle2d_basic() {
        let c = make_circle2d(DVec2::ZERO, 5.0).unwrap();
        assert!((c.radius - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_make_circle2d_3p() {
        let c = make_circle2d_3p(
            DVec2::new(1.0, 0.0),
            DVec2::new(0.0, 1.0),
            DVec2::new(-1.0, 0.0),
        )
        .unwrap();
        assert!((c.radius - 1.0).abs() < 1e-7);
        assert!((c.center - DVec2::ZERO).length() < 1e-7);
    }

    #[test]
    fn test_make_line2d_2p() {
        let l = make_line2d_2p(DVec2::ZERO, DVec2::new(1.0, 0.0)).unwrap();
        assert!((l.direction - DVec2::X).length() < 1e-12);
    }

    #[test]
    fn test_make_segment2d() {
        let seg = make_segment2d(DVec2::ZERO, DVec2::new(3.0, 0.0)).unwrap();
        if let Curve2d::Trimmed(tc) = &seg {
            assert!((tc.t_min - 0.0).abs() < 1e-12);
        } else {
            panic!("expected Trimmed curve");
        }
    }

    #[test]
    fn test_make_ellipse2d() {
        let e = make_ellipse2d(DVec2::ZERO, DVec2::X, 5.0, 3.0).unwrap();
        assert!((e.major_radius - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_make_parabola2d() {
        let p = make_parabola2d(DVec2::ZERO, DVec2::X, 2.0).unwrap();
        assert!((p.focal_param - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_make_arc_of_circle2d_params() {
        let circle = Circle2d::new(DVec2::ZERO, 5.0);
        let arc = make_arc_of_circle2d_params(&circle, 0.0, std::f64::consts::PI).unwrap();
        if let Curve2d::Trimmed(tc) = &arc {
            assert!((tc.t_min - 0.0).abs() < 1e-12);
            assert!((tc.t_max - std::f64::consts::PI).abs() < 1e-12);
        } else {
            panic!("expected Trimmed curve");
        }
    }
}
