//! 2D curve conversion utilities (Geom2dConvert).
//!
//! OCCT TKGeomBase Geom2dConvert package.
//!
//! Parallel to GeomConvert but for 2D parameter-space curves (Geom2d_*).
//! Provides BSpline↔Bezier conversion and curve composition.

use crate::geom::{BSplineCurve2, Curve2d};

pub mod bspline_curve;
pub mod bspline_curve_to_bezier_curve;

pub use bspline_curve::Geom2dBSplineCurve;
pub use bspline_curve_to_bezier_curve::BSplineCurveToBezierCurve;

/// Compose multiple 2D curves into a single BSpline.
///
/// OCCT: `Geom2dConvert::CompCurveToBSplineCurve`.
pub fn compose_curves_to_bspline(curves: &[Curve2d]) -> Option<BSplineCurve2> {
    use crate::geom::Curve2dEval;
    if curves.is_empty() {
        return None;
    }
    if curves.len() == 1 {
        return match &curves[0] {
            Curve2d::BSpline(b) => Some(b.clone()),
            _ => None,
        };
    }
    // Sample and interpolate
    let mut all_pts = Vec::new();
    for c in curves {
        let dom = c.default_domain();
        if !dom[0].is_finite() || !dom[1].is_finite() {
            return None;
        }
        let n = 8.max(2);
        for i in 0..n {
            let t = dom[0] + (dom[1] - dom[0]) * (i as f64) / ((n - 1) as f64);
            all_pts.push(c.point_at(t));
        }
    }
    if all_pts.len() < 2 {
        return None;
    }
    Some(BSplineCurve2::approximate(&all_pts))
}

/// Approximate a 2D curve to BSpline within tolerance.
///
/// OCCT: `Geom2dConvert::ApproxCurve`.
pub fn approx_curve_to_bspline(curve: &Curve2d, tol: f64) -> Option<BSplineCurve2> {
    use crate::geom::Curve2dEval;
    let dom = curve.default_domain();
    if !dom[0].is_finite() || !dom[1].is_finite() {
        return None;
    }
    let range = dom[1] - dom[0];
    let n = (range.abs() / (tol * 10.0)).ceil() as usize;
    let n = n.clamp(8, 200);

    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = dom[0] + range * (i as f64) / ((n - 1) as f64);
        pts.push(curve.point_at(t));
    }
    if pts.len() < 2 {
        return None;
    }
    Some(BSplineCurve2::approximate(&pts))
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;
    use glam::DVec2;

    #[test]
    fn test_approx_curve() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let bs = approx_curve_to_bspline(&c, 0.1);
        assert!(bs.is_some());
    }
}
