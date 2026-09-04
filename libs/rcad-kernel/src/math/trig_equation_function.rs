// OCCT math_TrigonometricEquationFunction (math_TrigonometricEquationFunction.hxx)
// 1:1 Rust translation.
//
// This is function, which corresponds trigonometric equation
// a*cos(x)*cos(x) + 2*b*cos(x)*sin(x) + c*cos(x) + d*sin(x) + e = 0
// See class math_TrigonometricFunctionRoots.

use super::root::FunctionWithDerivative;

/// OCCT math_TrigonometricEquationFunction.
#[derive(Debug, Clone, Copy)]
pub struct TrigonometricEquationFunction {
    my_aa: f64,
    my_bb: f64,
    my_cc: f64,
    my_dd: f64,
    my_ee: f64,
}

impl TrigonometricEquationFunction {
    /// OCCT math_TrigonometricEquationFunction(A, B, C, D, E).
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        TrigonometricEquationFunction {
            my_aa: a,
            my_bb: b,
            my_cc: c,
            my_dd: d,
            my_ee: e,
        }
    }
}

impl super::root::FunctionValue for TrigonometricEquationFunction {
    /// OCCT Value(X, F) (hxx L47-53).
    fn value(&mut self, x: f64) -> Option<f64> {
        let cn = x.cos();
        let sn = x.sin();
        // F = AA*CN*CN + 2*BB*CN*SN + CC*CN + DD*SN + EE (expanded form kept).
        let f = cn * (self.my_aa * cn + (self.my_bb + self.my_bb) * sn + self.my_cc)
            + self.my_dd * sn
            + self.my_ee;
        Some(f)
    }
}

impl FunctionWithDerivative for TrigonometricEquationFunction {
    /// OCCT Derivative(X, D) (hxx L55-63).
    fn derivative(&mut self, x: f64) -> Option<f64> {
        let cn = x.cos();
        let sn = x.sin();
        // D = -2*AA*CN*SN + 2*BB*(CN*CN - SN*SN) - CC*SN + DD*CN;
        // the `D += D;` below is the x2 expansion form — kept as written.
        let mut d = -self.my_aa * cn * sn + self.my_bb * (cn * cn - sn * sn);
        d += d;
        d += -self.my_cc * sn + self.my_dd * cn;
        Some(d)
    }

    /// OCCT Values(X, F, D) (hxx L65-78).
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let cn = x.cos();
        let sn = x.sin();
        let aacn = self.my_aa * cn;
        let bbsn = self.my_bb * sn;

        let f = aacn * cn + bbsn * (cn + cn) + self.my_cc * cn + self.my_dd * sn + self.my_ee;
        let mut d = -aacn * sn + self.my_bb * (cn * cn - sn * sn);
        d += d;
        d += -self.my_cc * sn + self.my_dd * cn;
        Some((f, d))
    }
}
