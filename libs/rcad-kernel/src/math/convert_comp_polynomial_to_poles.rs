//! OCCT Convert_CompPolynomialToPoles (TKMath Convert package).
//!
//! 1:1 translation of `Convert_CompPolynomialToPoles.cxx` — the 8-array
//! constructor (L89-138) and `Perform` (L178-261) used by
//! AdvApprox_ApproxAFunction: converts a piecewise polynomial (canonical
//! coefficients per segment) into a BSpline by interpolation at Schoenberg
//! points.

use crate::math::bspl_lib::{build_schoenberg_points, interpolate, knot_sequence};
use crate::math::plib::no_derivative_eval_polynomial_flat;

/// OCCT Convert_CompPolynomialToPoles (variable-continuity constructor,
/// L89-138).  Arrays are flat with OCCT 1-based index layouts:
/// - `continuity[1..=num_curves]` (entries 2..=num_curves are consumed)
/// - `num_coeff_per_curve[lower .. lower+num_curves-1]`
/// - `polynomial_intervals` row-major with 2 columns
/// - `true_intervals[1..=num_curves+1]`
#[allow(clippy::too_many_arguments)]
pub struct ConvertCompPolynomialToPoles {
    degree: usize,
    done: bool,
    knots: Vec<f64>,
    mults: Vec<i32>,
    /// Row-major (num_poles x dimension).
    poles: Vec<f64>,
    nb_poles: usize,
}

impl ConvertCompPolynomialToPoles {
    /// OCCT constructor #2 (Convert_CompPolynomialToPoles.cxx L89-138).
    #[allow(clippy::too_many_arguments)]
    pub fn from_arrays(
        num_curves: usize,
        dimension: usize,
        max_degree: usize,
        continuity: &[i32],
        num_coeff_per_curve: &[i32],
        coefficients: &[f64],
        polynomial_intervals: &[f64],
        true_intervals: &[f64],
    ) -> Self {
        assert!(
            num_curves > 0 && max_degree > 0 && dimension > 0 && polynomial_intervals.len() >= 2,
            "Convert_CompPolynomialToPoles:bad arguments"
        );

        let mut degree = 0usize;
        for ii in 0..num_curves {
            degree = degree.max((num_coeff_per_curve[ii] - 1) as usize);
        }

        let mut knots = vec![0.0f64; num_curves + 1];
        for ii in 0..=num_curves {
            knots[ii] = true_intervals[ii];
        }

        let mut mults = vec![0i32; num_curves + 1];
        for ii in 1..num_curves {
            // OCCT: if ((Continuity(ii) > myDegree) && (NumCurves > 1)) throw.
            assert!(
                !(continuity[ii] as usize > degree && num_curves > 1),
                "Convert_CompPolynomialToPoles:Continuity is too great"
            );
            mults[ii] = degree as i32 - continuity[ii];
        }
        mults[0] = degree as i32 + 1;
        mults[num_curves] = degree as i32 + 1;

        let mut this = ConvertCompPolynomialToPoles {
            degree,
            done: false,
            knots,
            mults,
            poles: Vec::new(),
            nb_poles: 0,
        };
        this.perform(
            num_curves,
            max_degree,
            dimension,
            num_coeff_per_curve,
            coefficients,
            polynomial_intervals,
            true_intervals,
        );
        this
    }

    /// OCCT Convert_CompPolynomialToPoles::Perform (L178-261).
    #[allow(clippy::too_many_arguments)]
    fn perform(
        &mut self,
        num_curves: usize,
        max_degree: usize,
        dimension: usize,
        num_coeff_per_curve: &[i32],
        coefficients: &[f64],
        polynomial_intervals: &[f64],
        true_intervals: &[f64],
    ) {
        let degree_i = self.degree as i32;
        let mut num_flat_knots = (2 * degree_i + 2) as usize;
        for ii in 1..(self.mults.len() - 1) {
            num_flat_knots += self.mults[ii] as usize;
        }
        let num_poles = num_flat_knots - self.degree - 1;

        let mut flat_knots = vec![0.0f64; num_flat_knots];
        knot_sequence(&self.knots, &self.mults, self.degree, false, &mut flat_knots);

        let mut parameters = vec![0.0f64; num_poles];
        build_schoenberg_points(self.degree, &flat_knots, &mut parameters);

        self.poles = vec![0.0f64; num_poles * dimension];
        self.nb_poles = num_poles;

        let mut contact_array = vec![0i32; num_poles];

        let mut index = 2usize; // 1-based
        let mut tindex = 1usize; // TrueIntervals.Lower() + 1
        let mut pindex = 0usize; // PolynomialIntervals.LowerRow() (flat row base)

        for ii in 0..num_poles {
            contact_array[ii] = 0;
            while parameters[ii] >= true_intervals[tindex] && index <= num_curves {
                index += 1;
                tindex += 1;
                pindex += 1;
            }
            //
            // normalized value so that it fits the original intervals for
            // the polynomial definition of the curves
            //
            let mut normalized_value = parameters[ii] - true_intervals[tindex - 1];
            normalized_value /= true_intervals[tindex] - true_intervals[tindex - 1];
            normalized_value = (1.0 - normalized_value) * polynomial_intervals[pindex * 2]
                + normalized_value * polynomial_intervals[pindex * 2 + 1];
            let coeff_index =
                (index - 2) * dimension * (max_degree.max(self.degree) + 1);

            let coefficient_array = &coefficients[coeff_index..];
            let deg = (num_coeff_per_curve[index - 2] - 1) as i32;

            // PLib::NoDerivativeEvalPolynomial(normalized_value, Deg,
            // Dimension, Deg*Dimension, coefficient_array[0],
            // poles_array[poles_index]).
            let out = &mut self.poles[ii * dimension..(ii + 1) * dimension];
            no_derivative_eval_polynomial_flat(
                normalized_value,
                deg,
                dimension as i32,
                deg * dimension as i32,
                coefficient_array,
                out,
            );
        }
        //
        // interpolation at schoenberg points should yield the desired result
        //
        let inversion_problem = interpolate(
            self.degree,
            &flat_knots,
            &parameters,
            &contact_array,
            dimension,
            &mut self.poles,
        );
        assert!(
            inversion_problem == 0,
            "Convert_CompPolynomialToPoles:inversion_problem"
        );
        self.done = true;
    }

    /// OCCT Convert_CompPolynomialToPoles::IsDone.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Convert_CompPolynomialToPoles::NbPoles.
    pub fn nb_poles(&self) -> usize {
        if self.done {
            return self.nb_poles;
        }
        0
    }

    /// OCCT Convert_CompPolynomialToPoles::Poles — row-major flat
    /// (num_poles x dimension).
    pub fn poles(&self) -> &[f64] {
        assert!(self.done, "Convert_CompPolynomialToPoles::Poles");
        &self.poles
    }

    /// OCCT Convert_CompPolynomialToPoles::Knots.
    pub fn knots_vec(&self) -> Vec<f64> {
        assert!(self.done, "Convert_CompPolynomialToPoles::Knots");
        self.knots.clone()
    }

    /// OCCT Convert_CompPolynomialToPoles::Multiplicities.
    pub fn multiplicities_vec(&self) -> Vec<i32> {
        assert!(self.done, "Convert_CompPolynomialToPoles::Multiplicities");
        self.mults.clone()
    }

    /// OCCT Convert_CompPolynomialToPoles::Degree.
    pub fn degree(&self) -> usize {
        if self.done {
            return self.degree;
        }
        0
    }
}
