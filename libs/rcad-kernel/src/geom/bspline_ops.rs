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
    eval_homogeneous, increase_degree_count_knots,
    move_point_and_tangent as bspl_move_point_and_tangent,
    increase_degree as bspl_increase_degree, rational_derivatives_inplace,
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

impl BSplineCurve3 {
    /// OCCT Geom_BSplineCurve::DN — the N-th derivative through the
    /// homogeneous evaluation (BSplCLib::Eval) and PLib::RationalDerivatives
    /// (in-place form).
    pub fn dn(&self, u: f64, n: usize) -> DVec3 {
        let dim = 3usize;
        let mut poles_flat = Vec::with_capacity(self.control_points.len() * dim);
        for p in &self.control_points {
            poles_flat.extend([p.x, p.y, p.z]);
        }
        let count = n + 1;
        let mut poles_res = vec![0.0f64; count * dim];
        let mut weights_res = vec![0.0f64; count];
        let mut extrap = [0i32; 2];
        // Clamp to the parameter range (the DN callers evaluate in range).
        let u_clamped = u
            .clamp(self.first_parameter(), self.last_parameter());
        eval_homogeneous(
            u_clamped,
            self.is_periodic,
            n as i32,
            &mut extrap,
            self.degree,
            &self.knots,
            dim,
            &poles_flat,
            &self.weights,
            &mut poles_res,
            &mut weights_res,
        );
        rational_derivatives_inplace(n as i32, dim, &mut poles_res, &mut weights_res);
        let off = n * dim;
        DVec3::new(poles_res[off], poles_res[off + 1], poles_res[off + 2])
    }
}

impl BSplineCurve3 {
    /// OCCT BSplCLib::Resolution(Poles, ArrayDimension, NumPoles, Weights,
    /// FlatKnots, Degree, Tolerance3D, UTolerance)
    /// (BSplCLib.cxx L4316-4820) — the parametric resolution of a BSpline.
    /// OCCT unrolls ArrayDimension 2/3/4 with identical arithmetic per
    /// coordinate k (ascending k); the generic loop below keeps the same
    /// accumulation order.  3D curve path (ArrayDimension = 3).
    pub fn bsplclib_resolution(&self, tolerance_3d: f64) -> f64 {
        // OCCT setup (BSplCLib.cxx L4317-4337).
        let degree = self.degree;
        let deg1 = degree + 1;
        let deg2 = (degree << 1) + 1;
        let fk = &self.knots; // FlatKnots (0-based flat sequence)
        let num_poles = fk.len() - deg1; // OCCT: FlatKnots.Length() - Deg1
        let num_poles_occt = self.control_points.len(); // OCCT NumPoles argument
        let max_derivative: f64;

        if self.is_rational() {
            // OCCT Weights branch (dim-3 form, BSplCLib.cxx L4440-4560).
            let wg = &self.weights;
            let mut min_weights = wg[0];
            for &w in wg.iter().take(num_poles_occt).skip(1) {
                if w < min_weights {
                    min_weights = w;
                }
            }

            let mut max_der = 0.0f64;
            for ii in 1..num_poles {
                let ii_index = ii % num_poles_occt;
                let ii_minus = (ii - 1) % num_poles_occt;
                let wg_ii_index = wg[ii_index];
                let wg_ii_minus = wg[ii_minus];
                let inverse = 1.0 / (fk[ii + degree] - fk[ii]);
                let lower = (ii - deg1).max(0);
                let upper = (deg2 + ii).min(num_poles);

                for jj in lower..upper {
                    let jj_index = jj % num_poles_occt;
                    let mut value = 0.0f64;
                    for kk in 0..3 {
                        let pa_jj = self.control_points[jj_index][kk];
                        let pa_ii = self.control_points[ii_index][kk];
                        let pa_mi = self.control_points[ii_minus][kk];
                        let mut factor =
                            ((pa_jj - pa_ii) * wg_ii_index) - ((pa_jj - pa_mi) * wg_ii_minus);
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                    }
                    value *= inverse;
                    if max_der < value {
                        max_der = value;
                    }
                }
            }
            max_derivative = max_der / min_weights;
        } else {
            // OCCT non-weighted branch (dim-3 form, L4664-4700).
            let mut max_der = 0.0f64;
            for ii in 1..num_poles {
                let ii_index = ii % num_poles_occt;
                let ii_minus = (ii - 1) % num_poles_occt;
                let inverse = 1.0 / (fk[ii + degree] - fk[ii]);
                let mut value = 0.0f64;
                for kk in 0..3 {
                    let mut factor =
                        self.control_points[ii_index][kk] - self.control_points[ii_minus][kk];
                    if factor < 0.0 {
                        factor = -factor;
                    }
                    value += factor;
                }
                value *= inverse;
                if max_der < value {
                    max_der = value;
                }
            }
            max_derivative = max_der;
        }

        let max_derivative = max_derivative * degree as f64;
        if max_derivative > f64::MIN_POSITIVE {
            // OCCT: UTolerance = Tolerance3D / max_derivative.
            tolerance_3d / max_derivative
        } else {
            // OCCT: UTolerance = Tolerance3D / RealSmall().
            tolerance_3d / f64::MIN_POSITIVE
        }
    }

    /// OCCT BSplCLib::Intervals(theKnots, theMults, theDegree, isPeriodic,
    /// theContinuity, theFirst, theLast, theTolerance, theIntervals)
    /// (BSplCLib.cxx L4824-4948).  Returns the number of intervals; when
    /// `intervals_out` is Some, fills it with the interval bounds.
    pub fn bsplclib_intervals(
        &self,
        continuity: i32,
        first: f64,
        last: f64,
        tolerance: f64,
        mut intervals_out: Option<&mut Vec<f64>>,
    ) -> usize {
        use crate::math::bspl_lib::{
            at, ati, first_uknot_index_mults, last_uknot_index_mults, locate_parameter_main,
        };

        let (the_knots, the_mults) = self.knots_mults();

        // Remove all knots with multiplicity less or equal than
        // (degree - continuity) except first and last (BSplCLib.cxx L4833-4845).
        let degree = self.degree as i32;
        let a_first_index = if self.is_periodic {
            1
        } else {
            first_uknot_index_mults(degree as usize, &the_mults)
        };
        let a_last_index = if self.is_periodic {
            the_knots.len() as i32
        } else {
            last_uknot_index_mults(degree as usize, &the_mults)
        };
        let mut a_new_knots: Vec<f64> = Vec::new();
        for an_index in a_first_index..=a_last_index {
            if ati(&the_mults, an_index) > (degree - continuity)
                || an_index == a_first_index
                || an_index == a_last_index
            {
                a_new_knots.push(at(&the_knots, an_index));
            }
        }
        let a_nb_new_knots = a_new_knots.len() as i32;

        // The range boundaries (BSplCLib.cxx L4848-4877).
        let mut a_cur_first = first;
        let mut a_cur_last = last;
        let mut a_period = 0.0f64;
        let mut a_first_period = 0i32;
        let mut a_last_period = 0i32;
        if self.is_periodic {
            let a_lower = the_knots[0];
            let an_upper = the_knots[the_knots.len() - 1];
            a_period = an_upper - a_lower;

            while a_cur_first < a_lower {
                a_cur_first += a_period;
                a_first_period -= 1;
            }
            while a_cur_last < a_lower {
                a_cur_last += a_period;
                a_last_period -= 1;
            }
            while a_cur_first >= an_upper {
                a_cur_first -= a_period;
                a_first_period += 1;
            }
            while a_cur_last >= an_upper {
                a_cur_last -= a_period;
                a_last_period += 1;
            }
        }

        // Locate the left and nearest knot for boundaries (L4879-4899) — the
        // LocateParameter variant without multiplicities.
        let mut an_index1 = 0i32;
        let mut an_index2 = 0i32;
        let mut a_dummy_double = 0.0f64;
        locate_parameter_main(
            &a_new_knots,
            a_cur_first,
            false,
            1,
            a_nb_new_knots,
            &mut an_index1,
            &mut a_dummy_double,
            0.0,
            1.0,
        );
        locate_parameter_main(
            &a_new_knots,
            a_cur_last,
            false,
            1,
            a_nb_new_knots,
            &mut an_index2,
            &mut a_dummy_double,
            0.0,
            1.0,
        );

        // The case when the beginning of the range coincides with the next knot.
        if an_index1 < a_nb_new_knots
            && (a_new_knots[an_index1 as usize] - a_cur_first).abs() < tolerance
        {
            an_index1 += 1;
        }
        // The case when the ending of the range coincides with the current knot.
        if a_nb_new_knots > 0
            && (a_new_knots[(an_index2 - 1) as usize] - a_cur_last).abs() < tolerance
        {
            an_index2 -= 1;
        }
        let a_nb_intervals = (an_index2 - an_index1 + 1
            + (a_last_period - a_first_period) * (a_nb_new_knots - 1)) as usize;

        // Fill the interval array (BSplCLib.cxx L4922-4945).
        if let Some(out) = intervals_out.as_deref_mut() {
            out.clear();
            if self.is_periodic && a_last_period != a_first_period {
                // Part from the beginning of range to the end of the first period.
                let mut i = an_index1;
                while i < a_nb_new_knots {
                    out.push(a_new_knots[i as usize] + a_first_period as f64 * a_period);
                    i += 1;
                }
                // Full periods.
                let mut a_period_num = a_first_period + 1;
                while a_period_num < a_last_period {
                    let mut i = 1;
                    while i < a_nb_new_knots {
                        out.push(a_new_knots[i as usize] + a_period_num as f64 * a_period);
                        i += 1;
                    }
                    a_period_num += 1;
                }
                // Part from the beginning of the last period to the end of range.
                let mut i = 1;
                while i <= an_index2 {
                    out.push(a_new_knots[i as usize] + a_last_period as f64 * a_period);
                    i += 1;
                }
            } else {
                let mut i = an_index1;
                while i <= an_index2 {
                    out.push(a_new_knots[i as usize] + a_first_period as f64 * a_period);
                    i += 1;
                }
            }
            // Update the first position and write the ending of the range.
            out[0] = first;
            out.push(last);
        }

        a_nb_intervals
    }
}
