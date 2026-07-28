//! Hermite interpolation (Hermit).
//!
//! OCCT TKGeomBase Hermit package.
//! Constructs a curve from position+derivative constraints at endpoints.
//!
//! P(t) = (2t³-3t²+1)*P0 + (t³-2t²+t)*V0 + (-2t³+3t²)*P1 + (t³-t²)*V1
//! where P0,P1 = endpoints, V0,V1 = tangent vectors.

use glam::DVec3;
use crate::geom::{BSplineCurve3, Curve3};

const TOL: f64 = 1e-12;

/// Build a cubic Hermite curve from endpoints and tangents.
///
/// Returns a degree-3 BSpline curve.
///
/// OCCT: `Hermit::Hermit(P0, V0, P1, V1)`.
pub fn hermit_curve(p0: DVec3, v0: DVec3, p1: DVec3, v1: DVec3) -> BSplineCurve3 {
    // Hermite control points expressed as cubic Bezier:
    // B0 = P0
    // B1 = P0 + V0/3
    // B2 = P1 - V1/3
    // B3 = P1
    let b0 = p0;
    let b1 = p0 + v0 / 3.0;
    let b2 = p1 - v1 / 3.0;
    let b3 = p1;

    BSplineCurve3 {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![b0, b1, b2, b3],
        weights: vec![],
    }
}

/// Evaluates a cubic Hermite curve at parameter t in [0, 1].
pub fn hermit_eval(p0: DVec3, v0: DVec3, p1: DVec3, v1: DVec3, t: f64) -> DVec3 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0; // P0 weight
    let h10 = t3 - 2.0 * t2 + t;          // V0 weight
    let h01 = -2.0 * t3 + 3.0 * t2;       // P1 weight
    let h11 = t3 - t2;                    // V1 weight
    h00 * p0 + h10 * v0 + h01 * p1 + h11 * v1
}

/// Evaluates the derivative of a cubic Hermite curve at parameter t in [0, 1].
pub fn hermit_deriv(p0: DVec3, v0: DVec3, p1: DVec3, v1: DVec3, t: f64) -> DVec3 {
    let t2 = t * t;
    let dh00 = 6.0 * t2 - 6.0 * t;
    let dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    let dh01 = -6.0 * t2 + 6.0 * t;
    let dh11 = 3.0 * t2 - 2.0 * t;
    dh00 * p0 + dh10 * v0 + dh01 * p1 + dh11 * v1
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_hermit_endpoints() {
        let bs = hermit_curve(DVec3::ZERO, DVec3::X, DVec3::new(1.0, 1.0, 0.0), DVec3::X);
        let p0 = bs.point_at(0.0);
        let p1 = bs.point_at(1.0);
        assert!((p0 - DVec3::ZERO).length() < TOL);
        assert!((p1 - DVec3::new(1.0, 1.0, 0.0)).length() < TOL);
    }
}
