//! `Geom_BSplineCurve` operations (TKG3d Geom package) used by the helix
//! pipeline: `IncreaseDegree` (Geom_BSplineCurve.cxx L243-287),
//! `MovePointAndTangent` (L1123-1158), `Translated` (Geom_Curve) and
//! `Transform` (pole transformation), plus the (knots, mults) <-> flat knot
//! sequence conversions.
//!
//! rcad's `BSplineCurve3` stores the flat (expanded) knot vector; OCCT
//! stores (knots, mults) and maintains `myFlatKnots`.  The two are
//! bijective, so the conversions here preserve OCCT semantics exactly.

use super::BSplineCurve3;
use crate::math::bspl_lib::{
    increase_degree_count_knots, move_point_and_tangent as bspl_move_point_and_tangent,
    increase_degree as bspl_increase_degree,
};
use crate::math::gp::Trsf;
use glam::DVec3;

/// OCCT Geom_BSplineCurve::MaxDegree() == BSplCLib::MaxDegree() == 25.
pub const BSPLINE_MAX_DEGREE: usize = 25;

impl BSplineCurve3 {
    /// OCCT Geom_BSplineCurve::FirstParameter.
    pub fn first_parameter(&self) -> f64 {
        self.knots[self.degree]
    }

    /// OCCT Geom_BSplineCurve::LastParameter.
    pub fn last_parameter(&self) -> f64 {
        self.knots[self.knots.len() - self.degree - 1]
    }

    /// Split the flat knot vector into OCCT (knots, mults) form.
    pub fn knots_mults(&self) -> (Vec<f64>, Vec<i32>) {
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for (i, k) in self.knots.iter().enumerate() {
            if i > 0 && *k == knots.last().copied().unwrap_or(f64::NAN) {
                *mults.last_mut().unwrap() += 1;
            } else {
                knots.push(*k);
                mults.push(1);
            }
        }
        (knots, mults)
    }

    /// Build the flat knot vector from OCCT (knots, mults).
    pub fn from_knots_mults(
        degree: usize,
        knots: Vec<f64>,
        mults: Vec<i32>,
        control_points: Vec<DVec3>,
    ) -> Self {
        let mut flat = Vec::new();
        for (k, m) in knots.iter().zip(mults.iter()) {
            for _ in 0..*m {
                flat.push(*k);
            }
        }
        let weights = vec![1.0; control_points.len()];
        BSplineCurve3 {
            degree,
            knots: flat,
            control_points,
            weights,
            is_periodic: false,
        }
    }

    /// OCCT Geom_BSplineCurve::IncreaseDegree(Degree) (L243-287).
    /// Periodic curves are not supported here (never produced by the helix
    /// pipeline; OCCT handles them through the same routine with wrapped
    /// knots).
    pub fn increase_degree(&mut self, degree: usize) {
        if degree == self.degree {
            return;
        }
        assert!(
            degree > self.degree && degree <= BSPLINE_MAX_DEGREE,
            "BSpline curve: IncreaseDegree: bad degree value"
        );
        assert!(!self.is_periodic, "IncreaseDegree: periodic not supported");

        let (knots, mults) = self.knots_mults();
        let nbknots = increase_degree_count_knots(self.degree, degree, false, &mults);

        let poles_flat: Vec<f64> = self
            .control_points
            .iter()
            .flat_map(|p| [p.x, p.y, p.z])
            .collect();
        let step = degree - self.degree;
        // OCCT npoles sizing: myPoles.Length() + Step * (ToK2 - FromK1)
        // (non-periodic clamped: ToK2 - FromK1 == nb_knots - 1).
        let mut new_poles = vec![0.0f64; (self.control_points.len() + step * (mults.len() - 1)) * 3];
        let mut new_knots = vec![0.0f64; nbknots];
        let mut new_mults = vec![0i32; nbknots];
        bspl_increase_degree(
            self.degree,
            degree,
            false,
            3,
            &poles_flat,
            &knots,
            &mults,
            &mut new_poles,
            &mut new_knots,
            &mut new_mults,
        );

        self.degree = degree;
        let new_count = new_poles.len() / 3;
        self.control_points = (0..new_count)
            .map(|i| DVec3::new(new_poles[i * 3], new_poles[i * 3 + 1], new_poles[i * 3 + 2]))
            .collect();
        self.weights = vec![1.0; new_count];
        let mut flat = Vec::new();
        for (k, m) in new_knots.iter().zip(new_mults.iter()) {
            for _ in 0..*m {
                flat.push(*k);
            }
        }
        self.knots = flat;
    }

    /// The flat knots as a plain slice (helper mirroring `myFlatKnots`).
    /// OCCT Geom_BSplineCurve::MovePointAndTangent(U, P, Tangent, Tolerance,
    /// StartingCondition, EndingCondition, ErrorStatus) (L1123-1158).
    /// On success (error status 0) the curve poles are replaced.
    #[allow(clippy::too_many_arguments)]
    pub fn move_point_and_tangent(
        &mut self,
        u: f64,
        p: DVec3,
        tangent: DVec3,
        tolerance: f64,
        starting_condition: i32,
        ending_condition: i32,
    ) -> i32 {
        assert!(
            !self.is_periodic,
            "MovePointAndTangent: periodic curves need SetNotPeriodic (no consumer in the pipeline)"
        );
        let poles_flat: Vec<f64> = self
            .control_points
            .iter()
            .flat_map(|v| [v.x, v.y, v.z])
            .collect();

        // Geom_Curve::D1(U, P0, delta_derivative).
        let p0 = crate::math::bspl::de_boor(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            u,
        );
        let d1 = self.derivative_at(u);
        let delta = p - p0;
        let delta_derivative = tangent - d1;

        let mut new_poles = vec![0.0f64; poles_flat.len()];
        // OCCT passes Weights() — always present in OCCT 8.0 (unit weights
        // for non-rational curves), so the homogeneous branch is taken.
        let weights = vec![1.0f64; self.control_points.len()];
        let error_status = bspl_move_point_and_tangent(
            u,
            3,
            &[delta.x, delta.y, delta.z],
            &[delta_derivative.x, delta_derivative.y, delta_derivative.z],
            tolerance,
            self.degree,
            starting_condition,
            ending_condition,
            &poles_flat,
            Some(&weights),
            &self.knots,
            &mut new_poles,
        );
        if error_status == 0 {
            let count = new_poles.len() / 3;
            self.control_points = (0..count)
                .map(|i| DVec3::new(new_poles[i * 3], new_poles[i * 3 + 1], new_poles[i * 3 + 2]))
                .collect();
        }
        error_status
    }

    /// OCCT Geom_Curve::Translated(P1, P2) — returns a copy translated by
    /// the vector (P2 - P1).
    #[must_use]
    pub fn translated(&self, p1: DVec3, p2: DVec3) -> BSplineCurve3 {
        let delta = p2 - p1;
        let mut out = self.clone();
        for v in out.control_points.iter_mut() {
            *v += delta;
        }
        out
    }

    /// OCCT Geom_BSplineCurve::Transform(theT) — transforms the poles.
    pub fn transform_trsf(&mut self, trsf: &Trsf) {
        for v in self.control_points.iter_mut() {
            *v = trsf.apply(*v);
        }
    }
}

impl BSplineCurve3 {
    /// OCCT Geom_BSplineCurve::IsRational — architecture note: rcad stores a
    /// weight array with an implicit 1.0 for non-rational curves, so
    /// rationality is observed as "any weight != 1".
    pub fn is_rational(&self) -> bool {
        self.weights.iter().any(|&w| w != 1.0)
    }

    /// OCCT Geom_BSplineCurve::Reverse (Geom_BSplineCurve.cxx L496-516):
    /// BSplCLib::Reverse(myKnots) + BSplCLib::Reverse(myMults) +
    /// BSplCLib::Reverse(myPoles, last) + BSplCLib::Reverse(myWeights, last),
    /// then the flat knot vector is rebuilt.  On rcad's flat knot vector,
    /// reversing (knots, mults) and re-expanding equals order-reversing the
    /// flat array while reflecting each value (kfirst + klast - k) — the
    /// combined BSplCLib::Reverse(Knots) / Reverse(Mults) effect.
    pub fn reversed(&self) -> BSplineCurve3 {
        let kfirst = self.knots[0];
        let klast = self.knots[self.knots.len() - 1];
        let knots: Vec<f64> = self.knots.iter().rev().map(|k| kfirst + klast - k).collect();
        let control_points: Vec<DVec3> = self.control_points.iter().rev().copied().collect();
        let mut weights: Vec<f64> = self.weights.clone();
        if self.is_rational() {
            weights.reverse();
        }
        BSplineCurve3 {
            degree: self.degree,
            knots,
            control_points,
            weights,
            is_periodic: self.is_periodic,
        }
    }

    /// OCCT Geom_BSplineCurve::SetKnots — replaces the (knots, mults) arrays
    /// keeping the pole count; rcad rebuilds the flat knot vector.
    pub fn set_knots(&mut self, knots: &[f64], mults: &[i32]) {
        assert_eq!(knots.len(), mults.len(), "SetKnots: knots/mults length mismatch");
        let mut flat =
            Vec::with_capacity(knots.iter().zip(mults.iter()).map(|(_, m)| *m as usize).sum());
        for (k, m) in knots.iter().zip(mults.iter()) {
            for _ in 0..*m {
                flat.push(*k);
            }
        }
        self.knots = flat;
    }
}
