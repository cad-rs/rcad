//! OCCT GeomFill_SnglrFunc (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_SnglrFunc.hxx (L60-99 members) + GeomFill_SnglrFunc.cxx (whole
//! file L26-159).
//!
//! Architecture mapping: the OCCT class is an `Adaptor3d_Curve` whose
//! "value" at U is `(D1 x D2) * ratio` — an auxiliary curve locating the
//! singularities of the base curve.  In rcad the base `Adaptor3d_Curve` is
//! a `Curve3`; the derivative requests map through `CurveEval::derivatives`
//! / `derivatives2` and, for the 4th/5th order needed by EvalD2/EvalD3,
//! through the BSpline DN kernel (only BSpline bases appear in the Frenet
//! unit).

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


    /// OCCT EvalD0 (L93-97): P = (D1 x D2) * ratio.
    pub fn eval_d0(&self, u: f64) -> DVec3 {
        let (_, d1, d2, _) = self.d3_parts(u);
        d1.cross(d2) * self.ratio
    }

    /// OCCT EvalD1 (L99-105): P = (D1 x D2) * ratio, D1 = (D1 x D3).
    pub fn eval_d1(&self, u: f64) -> (DVec3, DVec3) {
        let (_, d1, d2, d3) = self.d3_parts(u);
        let dc = d1 * self.ratio;
        (dc.cross(d2), dc.cross(d3))
    }

    /// OCCT EvalD2 (L107-116):
    /// P = (D1 x D2) * ratio, D1 = (D1 x D3) * ratio,
    /// D2 = ((D2 x D3) + (D1 x D4)) * ratio.
    pub fn eval_d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        let (_, d1, d2, d3) = self.d3_parts(u);
        let d4 = self.dn(u, 4);
        (
            (d1.cross(d2)) * self.ratio,
            d1.cross(d3) * self.ratio,
            (d2.cross(d3) + d1.cross(d4)) * self.ratio,
        )
    }

    /// OCCT EvalD3 (L118-130) — with the 5th derivative.
    #[allow(dead_code)]
    pub fn eval_d3(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        let (_, d1, d2, d3) = self.d3_parts(u);
        let d4 = self.dn(u, 4);
        let d5 = self.dn(u, 5);
        (
            (d1.cross(d2)) * self.ratio,
            d1.cross(d3) * self.ratio,
            (d2.cross(d3) + d1.cross(d4)) * self.ratio,
            (d1.cross(d5) + d2.cross(d4) * 2.0) * self.ratio,
        )
    }

    /// OCCT DN (L132-152) — orders 1..3 only (higher orders raise).
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

    /// The base curve's D3 evaluation: Adaptor3d_Curve::EvalD3 for the rcad
    /// Curve3.  BSpline bases go through the homogeneous DN kernel; the
    /// analytic types have closed-form evaluations only up to D2 in rcad,
    /// and their singularity search (GeomFill_Frenet::Init) stays
    /// anchor-out-of-scope.
    fn d3_parts(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        match &self.my_hcurve {
            Curve3::BSpline(bs) => {
                let p = bs.point_at(u);
                let d1 = bs.dn(u, 1);
                let d2 = bs.dn(u, 2);
                let d3 = bs.dn(u, 3);
                (p, d1, d2, d3)
            }
            _ => unimplemented!(
                "GeomFill_SnglrFunc D3 for non-BSpline bases is anchor-out-of-scope"
            ),
        }
    }
}
