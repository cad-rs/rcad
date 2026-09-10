//! OCCT Bisector_FunctionInter — F(u) = |PC(u)-PBis1(u)| - |PC(u)-PBis2(u)|,
//! 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_FunctionInter.hxx (L30-58) / .cxx (L29-107).

use std::sync::Arc;

use rcad_kernel::math::root::{FunctionValue, FunctionWithDerivative};

use super::bisector_curve::{BisectorCurve, GP_RESOLUTION};

/// OCCT Precision::Infinite().
const PRECISION_INFINITE: f64 = 2.0e100;

/// OCCT Bisector_FunctionInter (Bisector_FunctionInter.hxx L30-58).
pub struct FunctionInter {
    /// OCCT curve.
    curve: Option<Arc<dyn BisectorCurve>>,
    /// OCCT bisector1.
    bisector1: Option<Arc<dyn BisectorCurve>>,
    /// OCCT bisector2.
    bisector2: Option<Arc<dyn BisectorCurve>>,
}

impl Default for FunctionInter {
    /// OCCT Bisector_FunctionInter() (L29).
    fn default() -> Self {
        FunctionInter { curve: None, bisector1: None, bisector2: None }
    }
}

impl FunctionInter {
    /// OCCT Bisector_FunctionInter() (L29).
    pub fn new() -> Self {
        FunctionInter::default()
    }

    /// OCCT Bisector_FunctionInter(C, B1, B2) (L33-40).
    pub fn with_curves(
        c: Arc<dyn BisectorCurve>,
        b1: Arc<dyn BisectorCurve>,
        b2: Arc<dyn BisectorCurve>,
    ) -> Self {
        FunctionInter { curve: Some(c), bisector1: Some(b1), bisector2: Some(b2) }
    }

    /// OCCT Perform(C, B1, B2) (L44-51).
    pub fn perform(
        &mut self,
        c: Arc<dyn BisectorCurve>,
        b1: Arc<dyn BisectorCurve>,
        b2: Arc<dyn BisectorCurve>,
    ) {
        self.curve = Some(c);
        self.bisector1 = Some(b1);
        self.bisector2 = Some(b2);
    }
}

impl FunctionValue for FunctionInter {
    /// OCCT Value(X, F) (L55-64).
    fn value(&mut self, x: f64) -> Option<f64> {
        let pc = self.curve.as_ref().unwrap().value(x);
        let pb1 = self.bisector1.as_ref().unwrap().value(x);
        let pb2 = self.bisector2.as_ref().unwrap().value(x);

        let f = pc.distance(pb1) - pc.distance(pb2);

        Some(f)
    }
}

impl FunctionWithDerivative for FunctionInter {
    /// OCCT Derivative(X, D) (L68-72) — delegates to Values.
    fn derivative(&mut self, x: f64) -> Option<f64> {
        self.values(x).map(|(_f, d)| d)
    }

    /// OCCT Values(X, F, D) (L76-107).
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let (pc, tc) = self.curve.as_ref().unwrap().d1(x);
        let (pb1, tb1) = self.bisector1.as_ref().unwrap().d1(x);
        let (pb2, tb2) = self.bisector2.as_ref().unwrap().d1(x);

        let f1 = pc.distance(pb1);
        let f2 = pc.distance(pb2);
        let f = f1 - f2;

        let df1 = if f1.abs() < GP_RESOLUTION {
            PRECISION_INFINITE
        } else {
            ((pc.x - pb1.x) * (tc.x - tb1.x) + (pc.y - pb1.y) * (tc.y - tb1.y)) / f1
        };
        let df2 = if f2.abs() < GP_RESOLUTION {
            PRECISION_INFINITE
        } else {
            ((pc.x - pb2.x) * (tc.x - tb2.x) + (pc.y - pb2.y) * (tc.y - tb2.y)) / f2
        };
        let d = df1 - df2;

        Some((f, d))
    }
}
