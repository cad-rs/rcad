//! OCCT ChFi3d_SearchSing (TKFillet/ChFi3d) — 1:1 port of
//! ChFi3d_SearchSing.hxx (L28-59) + ChFi3d_SearchSing.cxx (whole file
//! L27-77).  Searches singularities on fillet:
//! F(t) = (C1(t) - C2(t)).(C1'(t) - C2'(t)).
//! Consumed by ChFi3d_FilBuilder::SplitSurf (ChFi3d_FilBuilder.cxx L2293).
//!
//! Architecture mapping: `class ChFi3d_SearchSing : public
//! math_FunctionWithDerivative` is expressed by implementing the
//! [`FunctionWithDerivative`] trait; the `occ::handle<Geom_Curve>` members
//! map to `&Curve3` references (the rcad Curve3 enum is the Geom_Curve
//! carrier); `myC1->D1/D2` map to the CurveEval accessors.

use rcad_kernel::geom::{Curve3, CurveEval as _};
use rcad_kernel::math::root::{FunctionValue, FunctionWithDerivative};

/// OCCT ChFi3d_SearchSing — searches singularities on fillet
/// (ChFi3d_SearchSing.hxx L36-38).
pub struct ChFi3dSearchSing<'a> {
    /// OCCT: occ::handle<Geom_Curve> myC1.
    my_c1: &'a Curve3,
    /// OCCT: occ::handle<Geom_Curve> myC2.
    my_c2: &'a Curve3,
}

impl<'a> ChFi3dSearchSing<'a> {
    /// OCCT ChFi3d_SearchSing(C1, C2) (ChFi3d_SearchSing.cxx L29-33).
    pub fn new(c1: &'a Curve3, c2: &'a Curve3) -> Self {
        ChFi3dSearchSing { my_c1: c1, my_c2: c2 }
    }
}

impl<'a> FunctionValue for ChFi3dSearchSing<'a> {
    /// OCCT Value(X, F) (ChFi3d_SearchSing.cxx L35-43) —
    /// gp_Vec V(P1, P2); F = V * (V2 - V1).
    fn value(&mut self, x: f64) -> Option<f64> {
        let p1 = self.my_c1.point_at(x);
        let v1 = self.my_c1.derivative_at(x);
        let p2 = self.my_c2.point_at(x);
        let v2 = self.my_c2.derivative_at(x);
        let v = p2 - p1;
        Some(v.dot(v2 - v1))
    }
}

impl<'a> FunctionWithDerivative for ChFi3dSearchSing<'a> {
    /// OCCT Derivative(X, D) (ChFi3d_SearchSing.cxx L45-56) —
    /// VPrim = V2 - V1; D = VPrim.SquareMagnitude() + (V * (W2 - W1)).
    fn derivative(&mut self, x: f64) -> Option<f64> {
        let p1 = self.my_c1.point_at(x);
        let v1 = self.my_c1.derivative_at(x);
        let w1 = self.my_c1.derivative2_at(x);
        let p2 = self.my_c2.point_at(x);
        let v2 = self.my_c2.derivative_at(x);
        let w2 = self.my_c2.derivative2_at(x);
        let v = p2 - p1;
        let vprim = v2 - v1;
        Some(vprim.length_squared() + v.dot(w2 - w1))
    }

    /// OCCT Values(X, F, D) (ChFi3d_SearchSing.cxx L58-69) —
    /// F = V * VPrim; D = VPrim.SquareMagnitude() + (V * (W2 - W1)).
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let p1 = self.my_c1.point_at(x);
        let v1 = self.my_c1.derivative_at(x);
        let w1 = self.my_c1.derivative2_at(x);
        let p2 = self.my_c2.point_at(x);
        let v2 = self.my_c2.derivative_at(x);
        let w2 = self.my_c2.derivative2_at(x);
        let v = p2 - p1;
        let vprim = v2 - v1;
        Some((
            v.dot(vprim),
            vprim.length_squared() + v.dot(w2 - w1),
        ))
    }
}
