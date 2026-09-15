//! OCCT GeomFill_SnglrFunc (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_SnglrFunc.hxx (L60-99 members) + GeomFill_SnglrFunc.cxx (whole
//! file L26-159).
//!
//! Architecture mapping: the OCCT class is an `Adaptor3d_Curve` whose
//! "value" at U is `(D1 x D2) * ratio` — an auxiliary curve locating the
//! singularities of the base curve.  In rcad the base `Adaptor3d_Curve` is
//! a `Curve3`; `myHCurve->EvalD3(theU)` maps through the base-curve D3
//! (the BSpline homogeneous DN kernel for BSpline bases) and
//! `myHCurve->EvalDN(theU, 4/5)` through the kernel `Geom_Curve::EvalDN`
//! union `curve_dn` (the ElCLib closed forms at any order for the
//! elementary kinds, exactly as `GeomAdaptor_Curve::EvalDN`).

use rcad_kernel::geom::curve_dn::curve_dn;
use rcad_kernel::geom::{Curve3, CurveEval};
use glam::DVec3;

use super::trihedron_law::{curve_first_parameter, curve_last_parameter};

/// OCCT GeomFill_SnglrFunc.
#[derive(Debug, Clone)]
pub struct SnglrFunc {
    my_hcurve: Curve3,
    ratio: f64,
}

impl SnglrFunc {
    /// OCCT GeomFill_SnglrFunc(HC) (L28-32).
    pub fn new(hcurve: Curve3) -> Self {
        SnglrFunc {
            my_hcurve: hcurve,
            ratio: 1.0,
        }
    }

    /// OCCT SetRatio (L38-41).
    pub fn set_ratio(&mut self, ratio: f64) {
        self.ratio = ratio;
    }

    /// OCCT FirstParameter (L43-46).
    pub fn first_parameter(&self) -> f64 {
        curve_first_parameter(&self.my_hcurve)
    }

    /// OCCT LastParameter (L48-51).
    pub fn last_parameter(&self) -> f64 {
        curve_last_parameter(&self.my_hcurve)
    }


    /// OCCT EvalD0 (L100-104): P = (D1 x D2) * ratio.
    pub fn eval_d0(&self, u: f64) -> DVec3 {
        let (_, d1, d2, _) = self.d3_parts(u);
        d1.cross(d2) * self.ratio
    }

    /// OCCT EvalD1 (L106-111): P = (D1 x D2) * ratio, D1 = (D1 x D3).
    pub fn eval_d1(&self, u: f64) -> (DVec3, DVec3) {
        let (_, d1, d2, d3) = self.d3_parts(u);
        let dc = d1 * self.ratio;
        (dc.cross(d2), dc.cross(d3))
    }

    /// OCCT EvalD2 (L113-120):
    /// P = (D1 x D2) * ratio, D1 = (D1 x D3) * ratio,
    /// D2 = ((D2 x D3) + (D1 x D4)) * ratio.
    pub fn eval_d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        let (_, d1, d2, d3) = self.d3_parts(u);
        // OCCT L116: aD4 = myHCurve->EvalDN(theU, 4) — the BASE curve's
        // 4th derivative (not the SnglrFunc wrapper's own DN).
        let d4 = self.base_dn(u, 4);
        (
            (d1.cross(d2)) * self.ratio,
            d1.cross(d3) * self.ratio,
            (d2.cross(d3) + d1.cross(d4)) * self.ratio,
        )
    }

    /// OCCT EvalD3 (L122-131) — with the 5th derivative.
    #[allow(dead_code)]
    pub fn eval_d3(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        let (_, d1, d2, d3) = self.d3_parts(u);
        let d4 = self.base_dn(u, 4);
        let d5 = self.base_dn(u, 5);
        (
            (d1.cross(d2)) * self.ratio,
            d1.cross(d3) * self.ratio,
            (d2.cross(d3) + d1.cross(d4)) * self.ratio,
            (d1.cross(d5) + d2.cross(d4) * 2.0) * self.ratio,
        )
    }

    /// OCCT DN (L133-152) — orders 1..3 only (higher orders raise).
    pub fn dn(&self, u: f64, n: usize) -> DVec3 {
        match n {
            1 => self.eval_d1(u).1,
            2 => self.eval_d2(u).2,
            3 => self.eval_d3(u).3,
            _ => panic!(
                "Exception: Derivative order is greater than 3. Cannot compute of derivative."
            ),
        }
    }

    /// OCCT myHCurve->EvalDN(theU, N) (cxx L116/L125-126) — the base
    /// curve's DN: the kernel `Geom_Curve::EvalDN` union over the rcad
    /// Curve3 ([`curve_dn`]; the ElCLib closed forms at any order for the
    /// elementary kinds — exactly the `GeomAdaptor_Curve::EvalDN` arms —
    /// the BSpline/Bezier DN kernel, the trimmed wrapper recursing on its
    /// basis).
    fn base_dn(&self, u: f64, n: usize) -> DVec3 {
        curve_dn(&self.my_hcurve, u, n as i32)
    }

    /// The base curve's D3 evaluation: OCCT `myHCurve->EvalD3(theU)`
    /// (cxx L102/L108/L115/L124 — the Adaptor3d_Curve D3 shared by
    /// EvalD0/D1/D2/D3).  The BSpline arm rides the homogeneous DN kernel
    /// (the rcad encoding of the OCCT BSpline D3, the same engine as
    /// `base_dn`); the other kinds take the kernel `Geom_Curve::EvalDN`
    /// union for the derivative slots (the ElCLib closed forms — the Line
    /// answers the null vector above D1, the Parabola the null vector
    /// above D2, as in the `GeomAdaptor_Curve::EvalD3` arms — and the
    /// trimmed wrapper recurses on its basis exactly as the
    /// `GeomAdaptor_Curve::load` unwrapping).
    fn d3_parts(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        match &self.my_hcurve {
            Curve3::BSpline(bs) => {
                let p = bs.point_at(u);
                let d1 = bs.dn(u, 1);
                let d2 = bs.dn(u, 2);
                let d3 = bs.dn(u, 3);
                (p, d1, d2, d3)
            }
            _ => {
                let p = CurveEval::point_at(&self.my_hcurve, u);
                (
                    p,
                    curve_dn(&self.my_hcurve, u, 1),
                    curve_dn(&self.my_hcurve, u, 2),
                    curve_dn(&self.my_hcurve, u, 3),
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle3, Line3};

    /// OCCT GeomFill_SnglrFunc::EvalD3 (L122-131) over a line basis: the
    /// base EvalD3 answers (P, D1, 0, 0) and EvalDN(4)/EvalDN(5) the null
    /// vector (ElCLib::LineDN), so the wrapper answers the null point and
    /// null derivative vectors.
    #[test]
    fn eval_d3_line_basis_delegates_to_the_base_derivatives() {
        let dir = DVec3::new(1.0, 2.0, -2.0).normalize();
        let func = SnglrFunc::new(Curve3::Line(Line3::new(DVec3::new(1.0, -1.0, 0.5), dir)));
        let (p, d1, d2, d3) = func.eval_d3(0.37);
        assert_eq!(p, DVec3::ZERO);
        assert_eq!(d1, DVec3::ZERO);
        assert_eq!(d2, DVec3::ZERO);
        assert_eq!(d3, DVec3::ZERO);
    }

    /// OCCT GeomFill_SnglrFunc::EvalD3 (L122-131) over a circle basis: the
    /// chain reads the base EvalD3 (ElCLib::CircleD3) and the base
    /// EvalDN(4)/EvalDN(5) (ElCLib::CircleDN — V4 = R(cos X + sin Y),
    /// V5 = V1), so the wrapper answers
    /// (D1 x D2, D1 x D3, (D2 x D3) + (D1 x D4), (D1 x D5) + 2 (D2 x D4))
    /// of those closed forms exactly.
    #[test]
    fn eval_d3_circle_basis_delegates_to_the_base_derivatives() {
        let radius = 2.5;
        let circle = Circle3::new(DVec3::new(1.0, 1.0, 1.0), DVec3::Z, radius);
        let func = SnglrFunc::new(Curve3::Circle(circle.clone()));
        let u = 0.7;
        let (p, d1, d2, d3) = func.eval_d3(u);

        // The base curve's closed-form derivatives (the ElCLib gp_Ax2
        // frame forms read through the circle payload).
        let (x, y) = (circle.x_dir, circle.y_dir);
        let v1 = radius * (-u.sin() * x + u.cos() * y);
        let v2 = radius * (-u.cos() * x - u.sin() * y);
        let v3 = radius * (u.sin() * x - u.cos() * y);
        let v4 = radius * (u.cos() * x + u.sin() * y); // N % 4 == 0 arm
        let v5 = radius * (-u.sin() * x + u.cos() * y); // (N - 1) % 4 == 0 arm

        let expected_p = v1.cross(v2);
        let expected_d1 = v1.cross(v3);
        let expected_d2 = v2.cross(v3) + v1.cross(v4);
        let expected_d3 = v1.cross(v5) + v2.cross(v4) * 2.0;

        assert!((p - expected_p).length() < 1e-12);
        assert!((d1 - expected_d1).length() < 1e-12);
        assert!((d2 - expected_d2).length() < 1e-12);
        assert!((d3 - expected_d3).length() < 1e-12);
    }

    /// OCCT GeomFill_SnglrFunc::EvalD2 (L113-120) over the same circle
    /// basis: the third slot is ((D2 x D3) + (D1 x D4)) * ratio of the base
    /// closed forms.
    #[test]
    fn eval_d2_circle_basis_delegates_to_the_base_derivatives() {
        let radius = 2.5;
        let circle = Circle3::new(DVec3::new(1.0, 1.0, 1.0), DVec3::Z, radius);
        let func = SnglrFunc::new(Curve3::Circle(circle.clone()));
        let u = -1.1;
        let (p, d1, d2) = func.eval_d2(u);

        let (x, y) = (circle.x_dir, circle.y_dir);
        let v1 = radius * (-u.sin() * x + u.cos() * y);
        let v2 = radius * (-u.cos() * x - u.sin() * y);
        let v3 = radius * (u.sin() * x - u.cos() * y);
        let v4 = radius * (u.cos() * x + u.sin() * y); // N % 4 == 0 arm

        assert!((p - v1.cross(v2)).length() < 1e-12);
        assert!((d1 - v1.cross(v3)).length() < 1e-12);
        assert!((d2 - (v2.cross(v3) + v1.cross(v4))).length() < 1e-12);
    }
}
