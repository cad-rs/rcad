//! Collect points on curves (CPnts).
//!
//! OCCT TKGeomBase CPnts package.
//! Sampled collection of points along a curve with uniform parameter spacing.
//!
//! Similar to GCPnts but simpler: just generates a point list without adaptive
//! refinement. Used for approximation and display.

use glam::DVec3;
use crate::geom::{Curve3, CurveEval};

/// Generate uniformly-spaced points along a curve.
///
/// OCCT: `CPnts_UniformDeflection` / `CPnts_AbscissaPoint`.
pub fn uniform_points(curve: &Curve3, n: usize) -> Vec<(f64, DVec3)> {
    let dom = curve.default_domain();
    if !dom[0].is_finite() || !dom[1].is_finite() || n < 2 {
        return vec![];
    }
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = dom[0] + (dom[1] - dom[0]) * (i as f64) / ((n - 1) as f64);
        pts.push((t, curve.point_at(t)));
    }
    pts
}

/// Generate points within a specific parameter range.
pub fn uniform_points_range(curve: &Curve3, t_min: f64, t_max: f64, n: usize) -> Vec<(f64, DVec3)> {
    if !t_min.is_finite() || !t_max.is_finite() || n < 2 {
        return vec![];
    }
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = t_min + (t_max - t_min) * (i as f64) / ((n - 1) as f64);
        pts.push((t, curve.point_at(t)));
    }
    pts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;
    #[test]
    fn test_uniform_points() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let pts = uniform_points(&line, 5);
        assert_eq!(pts.len(), 5);
    }
}
