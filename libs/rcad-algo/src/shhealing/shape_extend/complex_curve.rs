//! OCCT ShapeExtend_ComplexCurve (TKShHealing): `.hxx` L29-99, `.cxx`
//! L29-135 and `.lxx` L19-64 — a curve consisting of several segments with
//! the basic `Geom_Curve` interface implemented over them.
//!
//! Architecture mapping: the OCCT abstract base class maps to a Rust trait.
//! The five pure-virtuals (`NbCurves`, `Curve`, `LocateParameter`,
//! `LocalToGlobal`, `GetScaleFactor`) stay required; every concrete OCCT
//! method is provided with its 1:1 body.  The `myClosed` field (hxx L98)
//! maps to the `my_closed`/`set_my_closed` accessors (traits cannot hold
//! fields); `Curve(i)->Transform(T)` mutates the segment through the OCCT
//! handle, so the trait carries a `curve_mut` accessor as the bridge.
//!
//! Architecture mapping of the segments: `handle(Geom_Curve)` maps to the
//! rcad value type [`rcad_kernel::geom::Curve3`]; evaluation uses the
//! `CurveEval` trait (`point_at` = Value/D0, `derivative_at` = D1,
//! `derivative2_at` = D2, `derivative3_at` = D3, `default_domain` =
//! First/LastParameter).

use rcad_kernel::geom::{transform_curve, CurveEval, Curve3};
use rcad_kernel::math::GeomAbsShape;
use glam::DVec3;

/// OCCT Geom_Curve::ResD1 (Geom_Curve.hxx L62-66): point and first derivative.
#[derive(Debug, Clone, Copy)]
pub struct CurveResD1 {
    pub point: DVec3,
    pub d1: DVec3,
}

/// OCCT Geom_Curve::ResD2 (Geom_Curve.hxx L69-74): point and first two
/// derivatives.
#[derive(Debug, Clone, Copy)]
pub struct CurveResD2 {
    pub point: DVec3,
    pub d1: DVec3,
    pub d2: DVec3,
}

/// OCCT Geom_Curve::ResD3 (Geom_Curve.hxx L77-83): point and first three
/// derivatives.
#[derive(Debug, Clone, Copy)]
pub struct CurveResD3 {
    pub point: DVec3,
    pub d1: DVec3,
    pub d2: DVec3,
    pub d3: DVec3,
}

/// OCCT ShapeExtend_ComplexCurve (ShapeExtend_ComplexCurve.hxx L29-99).
pub trait ShapeExtendComplexCurve {
    /// OCCT NbCurves() (hxx L34): returns number of curves.  Pure virtual.
    fn nb_curves(&self) -> i32;

    /// OCCT Curve(index) (hxx L37): returns curve given by its index.
    /// Pure virtual.  The handle maps to an owned Curve3 clone.
    fn curve(&self, index: i32) -> Curve3;

    /// OCCT LocateParameter(U, UOut) (hxx L41): returns number of the curve
    /// for the given parameter U and local parameter UOut for the found
    /// curve.  Pure virtual; the C++ output parameter maps to the second
    /// tuple element.
    fn locate_parameter(&self, u: f64) -> (i32, f64);

    /// OCCT LocalToGlobal(index, Ulocal) (hxx L45): returns global parameter
    /// for the whole curve according to the segment and local parameter on
    /// it.  Pure virtual.
    fn local_to_global(&self, index: i32, ulocal: f64) -> f64;

    /// OCCT GetScaleFactor(ind) (hxx L84): returns scale factor for
    /// recomputing of deviatives.  Pure virtual.
    fn get_scale_factor(&self, ind: i32) -> f64;

    /// OCCT myClosed field (hxx L98) — trait field bridge (getter).
    fn my_closed(&self) -> bool;

    /// OCCT myClosed field (hxx L98) — trait field bridge (setter).
    fn set_my_closed(&mut self, closed: bool);

    /// Handle-mutation bridge for OCCT `Curve(i)->Transform(T)`: shared
    /// access to the i-th segment for in-place transformation.
    fn curve_mut(&mut self, index: i32) -> &mut Curve3;

    /// OCCT Transform(T) (cxx L36-42): applies transformation to each curve.
    fn transform(&mut self, t: &glam::DAffine3) {
        for i in 1..=self.nb_curves() {
            let transformed = transform_curve(self.curve_mut(i), t);
            *self.curve_mut(i) = transformed;
        }
    }

    /// OCCT EvalD0(U) (cxx L46-51): returns point at parameter U; finds the
    /// appropriate curve and local parameter on it.
    fn eval_d0(&self, u: f64) -> DVec3 {
        let (an_ind, u_out) = self.locate_parameter(u);
        self.curve(an_ind).point_at(u_out)
    }

    /// OCCT EvalD1(U) (cxx L55-62).
    fn eval_d1(&self, u: f64) -> CurveResD1 {
        let (an_ind, u_out) = self.locate_parameter(u);
        let c = self.curve(an_ind);
        let mut a_result = CurveResD1 {
            point: c.point_at(u_out),
            d1: c.derivative_at(u_out),
        };
        self.transform_dn(&mut a_result.d1, an_ind, 1);
        a_result
    }

    /// OCCT EvalD2(U) (cxx L66-74).
    fn eval_d2(&self, u: f64) -> CurveResD2 {
        let (an_ind, u_out) = self.locate_parameter(u);
        let c = self.curve(an_ind);
        let mut a_result = CurveResD2 {
            point: c.point_at(u_out),
            d1: c.derivative_at(u_out),
            d2: c.derivative2_at(u_out),
        };
        self.transform_dn(&mut a_result.d1, an_ind, 1);
        self.transform_dn(&mut a_result.d2, an_ind, 2);
        a_result
    }

    /// OCCT EvalD3(U) (cxx L78-87).
    fn eval_d3(&self, u: f64) -> CurveResD3 {
        let (an_ind, u_out) = self.locate_parameter(u);
        let c = self.curve(an_ind);
        let mut a_result = CurveResD3 {
            point: c.point_at(u_out),
            d1: c.derivative_at(u_out),
            d2: c.derivative2_at(u_out),
            d3: c.derivative3_at(u_out),
        };
        self.transform_dn(&mut a_result.d1, an_ind, 1);
        self.transform_dn(&mut a_result.d2, an_ind, 2);
        self.transform_dn(&mut a_result.d3, an_ind, 3);
        a_result
    }

    /// OCCT EvalDN(U, N) (cxx L91-101).
    ///
    /// Architecture note: rcad `CurveEval` exposes D1-D3; orders beyond 3 map
    /// to the zero vector before the OCCT `if (N)` TransformDN guard.
    fn eval_dn(&self, u: f64, n: i32) -> DVec3 {
        let (an_ind, u_out) = self.locate_parameter(u);
        let c = self.curve(an_ind);
        let mut a_result = match n {
            1 => c.derivative_at(u_out),
            2 => c.derivative2_at(u_out),
            3 => c.derivative3_at(u_out),
            _ => DVec3::ZERO,
        };
        if n != 0 {
            self.transform_dn(&mut a_result, an_ind, n);
        }
        a_result
    }

    /// OCCT CheckConnectivity(Preci) (cxx L105-124): checks geometrical
    /// connectivity of the curves, including closure (sets field myClosed).
    fn check_connectivity(&mut self, preci: f64) -> bool {
        let nb_c = self.nb_curves();
        let mut ok = true;
        for i in 1..nb_c {
            if i == 1 {
                let p_first = self.eval_d0(self.first_parameter());
                let p_last = self.eval_d0(self.last_parameter());
                self.set_my_closed(p_first.distance(p_last) <= preci);
            }
            let ci = self.curve(i);
            let ci1 = self.curve(i + 1);
            ok &= ci.point_at(ci.default_domain()[1]).distance(ci1.point_at(ci1.default_domain()[0])) <= preci;
        }
        ok
    }

    /// OCCT TransformDN(V, ind, N) (cxx L128-135): transforms the derivative
    /// according to its order.
    fn transform_dn(&self, v: &mut DVec3, ind: i32, n: i32) {
        let fact = self.get_scale_factor(ind);
        for _i in 1..=n {
            *v *= fact;
        }
    }

    /// OCCT ReversedParameter(U) (lxx L19-22): returns 1 - U.
    fn reversed_parameter(&self, u: f64) -> f64 {
        1.0 - u
    }

    /// OCCT FirstParameter() (lxx L26-29): returns 0.
    fn first_parameter(&self) -> f64 {
        0.0
    }

    /// OCCT LastParameter() (lxx L33-36): returns 1.
    fn last_parameter(&self) -> f64 {
        1.0
    }

    /// OCCT IsClosed() (lxx L40-43).
    fn is_closed(&self) -> bool {
        self.my_closed()
    }

    /// OCCT IsPeriodic() (lxx L47-50): returns False.
    fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT Continuity() (lxx L54-57): returns GeomAbs_C0.
    fn continuity(&self) -> GeomAbsShape {
        GeomAbsShape::C0
    }

    /// OCCT IsCN(N) (lxx L61-64): returns False if N > 0.
    fn is_cn(&self, n: i32) -> bool {
        n <= 0
    }
}
