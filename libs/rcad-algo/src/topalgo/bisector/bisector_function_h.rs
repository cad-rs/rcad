//! OCCT Bisector_FunctionH — H(v) = (T1.P2(v) - P1)*||T(v)|| -
//! (T(v).P2(v) - P1)*||T1||, 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_FunctionH.hxx (L34-56) / .cxx (L24-85).

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::math::root::{FunctionValue, FunctionWithDerivative};

use super::bisector_curve::BisectorCurve;

/// OCCT Bisector_FunctionH (Bisector_FunctionH.hxx L34-56).
pub struct FunctionH {
    /// OCCT curve2 (the owned handle).
    curve2: Arc<dyn BisectorCurve>,
    /// OCCT p1 (t1 kept unit, mirroring the constructor normalization).
    p1: DVec2,
    t1: DVec2,
}

impl FunctionH {
    /// OCCT Bisector_FunctionH(C2, P1, T1) (L24-32) — T1 is normalized.
    pub fn new(c2: Arc<dyn BisectorCurve>, p1: DVec2, t1: DVec2) -> Self {
        let t1 = t1.normalize_or_zero();
        FunctionH { curve2: c2, p1, t1 }
    }
}

impl FunctionValue for FunctionH {
    /// OCCT Value(X, F) (L39-52): F = (p1 - P2) . (||T2||*T1 - T2).
    fn value(&mut self, x: f64) -> Option<f64> {
        // point sur C2 / tangente a C2 en V.
        let (p2, t2) = self.curve2.d1(x);

        let norm_t2 = t2.length();
        let ax = norm_t2 * self.t1.x - t2.x;
        let ay = norm_t2 * self.t1.y - t2.y;

        let f = (self.p1.x - p2.x) * ax + (self.p1.y - p2.y) * ay;

        Some(f)
    }
}

impl FunctionWithDerivative for FunctionH {
    /// OCCT Derivative(X, D) (L56-60) — delegates to Values.
    fn derivative(&mut self, x: f64) -> Option<f64> {
        self.values(x).map(|(_f, d)| d)
    }

    /// OCCT Values(X, F, D) (L64-85).
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        // point sur C2 / tangente a C2 / derivee seconde a C2 en V.
        let (p2, t2, t2v) = self.curve2.d2(x);

        let norm_t2 = t2.length();
        let ax = norm_t2 * self.t1.x - t2.x;
        let ay = norm_t2 * self.t1.y - t2.y;

        let f = (self.p1.x - p2.x) * ax + (self.p1.y - p2.y) * ay;

        let scal = t2.dot(t2v) / norm_t2;
        let d_ax = scal * self.t1.x - t2v.x;
        let d_ay = scal * self.t1.y - t2v.y;

        let d = -t2.x * ax
            - t2.y * ay
            + (self.p1.x - p2.x) * d_ax
            + (self.p1.y - p2.y) * d_ay;

        Some((f, d))
    }
}
