//! Curve and surface local properties — curvature extrema and inflections (LProp).
//!
//! OCCT TKGeomBase LProp package.

use crate::geom::{Curve2d, Curve2dEval};

const TOL: f64 = 1e-7;

/// Result of a curvature analysis.
#[derive(Debug, Clone)]
pub struct LPropResult {
    pub inflection_params: Vec<f64>,
    pub curvature_extrema: Vec<f64>,
}

/// Find inflection points and curvature extrema on a 2D curve.
pub fn curve_analysis_2d(curve: &Curve2d) -> LPropResult {
    let dom = curve.default_domain();
    if !dom[0].is_finite() || !dom[1].is_finite() {
        return LPropResult { inflection_params: vec![], curvature_extrema: vec![] };
    }

    const N: usize = 101;
    let mut inflections = Vec::new();
    let mut extrema = Vec::new();
    let dt = (dom[1] - dom[0]) / (N as f64);

    let mut curvatures = Vec::with_capacity(N + 1);
    for i in 0..=N {
        let t = dom[0] + dt * (i as f64);
        curvatures.push(curve.curvature_at(t));
    }

    for i in 1..=N {
        let k_prev = curvatures[i - 1];
        let k = curvatures[i];
        // Sign change → inflection
        if k_prev * k < 0.0 {
            let frac = k_prev.abs() / (k_prev.abs() + k.abs());
            let infl = dom[0] + dt * ((i as f64) - 1.0 + frac);
            inflections.push(infl);
        }
        // Local extremum
        if i > 1 && i < N {
            let k_pp = curvatures[i - 2];
            let k_n = curvatures[i + 1];
            if (k_prev > k_pp && k_prev > k_n) || (k_prev < k_pp && k_prev < k_n) {
                extrema.push(dom[0] + dt * ((i - 1) as f64));
            }
        }
    }

    inflections.dedup_by(|a, b| (*a - *b).abs() < TOL);
    extrema.dedup_by(|a, b| (*a - *b).abs() < TOL);
    LPropResult { inflection_params: inflections, curvature_extrema: extrema }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;
    #[test]
    fn test_line_has_no_inflections() {
        let line = Curve2d::Line(Line2d::new(DVec2::ZERO, DVec2::X));
        let r = curve_analysis_2d(&line);
        assert!(r.inflection_params.is_empty());
    }
}
