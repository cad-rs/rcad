//! OCCT PLib_JacobiPolynomial (TKMath PLib package).
//!
//! 1:1 translation of `PLib_JacobiPolynomial.cxx` (L52-372: constructor,
//! Points, Weights, MaxValue, MaxError, ReduceDegree, AverageError,
//! ToCoefficients).
//!
//! Not translated (no consumer in the ported pipeline; the basis-evaluation
//! cache lives in OCCT `PLib_JacobiPolynomial_Coeffs.pxx`): D0/D1/D2/D3 and
//! the D0123 basis recursion (FORTRAN MPOJAC port).

use super::gauss_points::gauss_points;
use super::p_lib_jacobi_data::{
    MAX_VALUES_DB_C0, MAX_VALUES_DB_C1, MAX_VALUES_DB_C2, TRANS_MATRIX_C0, TRANS_MATRIX_C1,
    TRANS_MATRIX_C2, WEIGHTS_DB0_C0, WEIGHTS_DB0_C1, WEIGHTS_DB0_C2, WEIGHTS_DB_C0, WEIGHTS_DB_C1,
    WEIGHTS_DB_C2,
};
use super::GeomAbsShape;
use super::MatD;
use super::VecD;

/// OCCT PLib::NivConstr (PLib.cxx L2196+): GeomAbs_Shape -> integer level.
pub fn niv_constr(constraint_order: GeomAbsShape) -> usize {
    match constraint_order {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::C1 => 1,
        GeomAbsShape::C2 => 2,
        _ => panic!("PLib::NivConstr - invalid ConstraintOrder"),
    }
}

/// The possible values of NbGaussPoints (PLib_JacobiPolynomial.cxx L32-34).
const NB_GAUSS_POINTS: [usize; 9] = [8, 10, 15, 20, 25, 30, 40, 50, 61];
/// OCCT THE_INVALID_VALUE.
const THE_INVALID_VALUE: f64 = -999.0;
/// Maximum supported polynomial degree.
const THE_MAX_DEGREE: usize = 30;

/// OCCT PLib_JacobiPolynomial — Jacobi polynomial basis relative to an order
/// of constraint.  P(t) = R(t) + W(t) * Q(t), W(t) = (1-t^2)^(2*nivConstr+2).
#[derive(Debug, Clone)]
pub struct JacobiPolynomial {
    work_degree: usize,
    niv_constr: usize,
    degree: usize,
}

impl JacobiPolynomial {
    /// OCCT PLib_JacobiPolynomial(theWorkDegree, theConstraintOrder).
    pub fn new(work_degree: usize, constraint_order: GeomAbsShape) -> Self {
        let niv_constr = niv_constr(constraint_order);
        let degree = work_degree - 2 * (niv_constr + 1);
        assert!(
            work_degree >= 2 * (niv_constr + 1),
            "WorkDegree too small for given ConstraintOrder"
        );
        assert!(degree <= THE_MAX_DEGREE, "Invalid Degree");
        JacobiPolynomial {
            work_degree,
            niv_constr,
            degree,
        }
    }

    /// OCCT PLib_JacobiPolynomial::WorkDegree.
    pub fn work_degree(&self) -> usize {
        self.work_degree
    }

    /// OCCT PLib_JacobiPolynomial::NivConstr.
    pub fn niv_constr(&self) -> usize {
        self.niv_constr
    }

    /// OCCT PLib_JacobiPolynomial::Points — the positive Legendre roots by
    /// increasing value.  `tab_points` has bounds (0, nb_gauss_points/2).
    pub fn points(&self, nb_gauss_points: usize, tab_points: &mut [f64]) {
        let valid = NB_GAUSS_POINTS.contains(&nb_gauss_points) && nb_gauss_points > self.degree;
        assert!(valid, "Invalid NbGaussPoints");

        let mut decreasing = VecD::new(nb_gauss_points);
        gauss_points(nb_gauss_points, &mut decreasing);

        // theTabPoints consist of only positive increasing values.
        for i in 1..=nb_gauss_points / 2 {
            tab_points[i] = decreasing.get(nb_gauss_points / 2 - i + 1);
        }
        tab_points[0] = if nb_gauss_points % 2 == 1 {
            0.0
        } else {
            THE_INVALID_VALUE
        };
    }

    /// OCCT PLib_JacobiPolynomial::Weights — Gauss weights for the positive
    /// roots.  `tab_weights` is (0..nb_gauss_points/2) x (0..degree).
    pub fn weights(&self, nb_gauss_points: usize, tab_weights: &mut MatD) {
        let niv = self.niv_constr;
        // Zero-initialize entire output array.
        for i in 0..tab_weights.n_rows() {
            for j in 0..tab_weights.n_cols() {
                tab_weights.set(i + 1, j + 1, 0.0);
            }
        }

        let weights_db: &[f64] = match niv {
            0 => &WEIGHTS_DB_C0,
            1 => &WEIGHTS_DB_C1,
            _ => &WEIGHTS_DB_C2,
        };
        let min_degree = 2 * (niv + 1);

        // Calculate offset into the weights database (OCCT pointer arithmetic
        // converted to a flat running index).
        let mut db_index = 0usize;
        if nb_gauss_points > 8 {
            db_index += 8 * (8 - min_degree) / 2;
        }
        if nb_gauss_points > 10 {
            db_index += 10 * (10 - min_degree) / 2;
        }
        if nb_gauss_points > 15 {
            db_index += ((15 - 1) / 2) * (15 - min_degree);
        }
        if nb_gauss_points > 20 {
            db_index += 20 * (20 - min_degree) / 2;
        }
        if nb_gauss_points > 25 {
            db_index += ((25 - 1) / 2) * (25 - min_degree);
        }
        if nb_gauss_points > 30 {
            db_index += 30 * (30 - min_degree) / 2;
        }
        if nb_gauss_points > 40 {
            db_index += 40 * (40 - min_degree) / 2;
        }
        if nb_gauss_points > 50 {
            db_index += 50 * (50 - min_degree) / 2;
        }

        let half_points = nb_gauss_points / 2;
        // Copy TabWeightsDB into theTabWeights (row-major per OCCT read order:
        // for j in columns, for i in rows).
        for j in 0..=self.degree {
            for i in 1..=half_points {
                tab_weights.set(i + 1, j + 1, weights_db[db_index]);
                db_index += 1;
            }
        }

        if nb_gauss_points % 2 == 1 {
            // Odd — fill row 0 with special values.
            let weights_db0: &[f64] = match niv {
                0 => &WEIGHTS_DB0_C0,
                1 => &WEIGHTS_DB0_C1,
                _ => &WEIGHTS_DB0_C2,
            };
            let mut db0_index = 0usize;
            if nb_gauss_points > 15 {
                db0_index += (15 - 1 - min_degree) / 2 + 1;
            }
            if nb_gauss_points > 25 {
                db0_index += (25 - 1 - min_degree) / 2 + 1;
            }
            let mut j = 0usize;
            while j <= self.degree {
                tab_weights.set(1, j + 1, weights_db0[db0_index]);
                db0_index += 1;
                j += 2;
            }
        } else {
            // Even — row 0 = THE_INVALID_VALUE.
            for j in 0..=self.degree {
                tab_weights.set(1, j + 1, THE_INVALID_VALUE);
            }
        }
    }

    /// OCCT PLib_JacobiPolynomial::MaxValue — fills `tab_max` from the
    /// database for the constraint level.
    pub fn max_value(&self, tab_max: &mut [f64]) {
        let db: &[f64] = match self.niv_constr {
            0 => &MAX_VALUES_DB_C0,
            1 => &MAX_VALUES_DB_C1,
            _ => &MAX_VALUES_DB_C2,
        };
        for (i, v) in tab_max.iter_mut().enumerate() {
            *v = db[i];
        }
    }

    /// OCCT PLib_JacobiPolynomial::MaxError — maximum error of W(t)Q(t)
    /// obtained by missing the coefficients of jac_coeff from new_degree+1 to
    /// work_degree.  `jac_coeff` is the flat coefficient array laid out as
    /// [coeff_idx * dimension + dim].
    pub fn max_error(&self, dimension: usize, jac_coeff: &[f64], new_degree: usize) -> f64 {
        let mut tab_max = vec![0.0f64; self.degree + 2];
        self.max_value(&mut tab_max);

        let beg_idx = 2 * (self.niv_constr + 1);
        let cut_idx = beg_idx.max(new_degree + 1);

        let mut max_err_dim = vec![0.0f64; dimension];
        for dim_idx in 1..=dimension {
            for coeff_idx in cut_idx..=self.work_degree {
                let coeff_value = jac_coeff[coeff_idx * dimension + dim_idx - 1];
                let basis_max = tab_max[coeff_idx - beg_idx];
                max_err_dim[dim_idx - 1] += coeff_value.abs() * basis_max;
            }
        }

        // math_Vector::Norm.
        max_err_dim.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    /// OCCT PLib_JacobiPolynomial::ReduceDegree — computes
    /// new_degree <= max_degree so that max_error <= tol.
    pub fn reduce_degree(
        &self,
        dimension: usize,
        max_degree: usize,
        tol: f64,
        jac_coeff: &[f64],
        new_degree: &mut usize,
        max_error: &mut f64,
    ) {
        let idx = 2 * (self.niv_constr + 1) - 1;
        let cut_idx = idx + 1;

        let mut max_err_dim = vec![0.0f64; dimension];
        let mut tab_max = vec![0.0f64; self.degree + 2];
        self.max_value(&mut tab_max);

        *new_degree = idx;
        *max_error = 0.0;

        // Search for theNewDegree from high degree to low.
        for i in (cut_idx..=self.work_degree).rev() {
            for idim in 1..=dimension {
                max_err_dim[idim - 1] += jac_coeff[i * dimension + idim - 1].abs() * tab_max[i - cut_idx];
            }
            let error = max_err_dim.iter().map(|v| v * v).sum::<f64>().sqrt();
            if error > tol && i <= max_degree {
                *new_degree = i;
                break;
            }
            *max_error = error;
        }

        // Fallback: find last non-negligible coefficient.
        if *new_degree == idx {
            const EPS: f64 = 1.0e-9;
            *new_degree = 0;
            for i in (1..=idx).rev() {
                let mut bid = 0.0f64;
                for idim in 1..=dimension {
                    bid += jac_coeff[i * dimension + idim - 1].abs();
                }
                if bid > EPS {
                    *new_degree = i;
                    break;
                }
            }
        }
    }

    /// OCCT PLib_JacobiPolynomial::AverageError.
    pub fn average_error(&self, dimension: usize, jac_coeff: &[f64], new_degree: usize) -> f64 {
        let cut_idx = (2 * (self.niv_constr + 1) + 1).max(new_degree + 1);
        let mut average_err = 0.0f64;

        for idim in 1..=dimension {
            for i in cut_idx..=self.degree {
                let c = jac_coeff[i * dimension + idim - 1];
                average_err += c * c;
            }
        }

        (average_err / 2.0).sqrt()
    }

    /// OCCT PLib_JacobiPolynomial::ToCoefficients — converts
    /// P(t) = R(t) + W(t) Q(t) into the canonical base.
    /// `jac_coeff` / `coefficients` are flat arrays laid out as
    /// [degree_idx * dimension + dim].
    pub fn to_coefficients(
        &self,
        dimension: usize,
        degree: usize,
        jac_coeff: &[f64],
        coefficients: &mut [f64],
    ) {
        let max_m = THE_MAX_DEGREE + 1;
        let half_degree = degree / 2;
        let double_dim = 2 * dimension;
        let trans: &[[f64; 496]; 2] = match self.niv_constr {
            0 => &TRANS_MATRIX_C0,
            1 => &TRANS_MATRIX_C1,
            _ => &TRANS_MATRIX_C2,
        };
        // OCCT THE_TRANS_MATRIX[niv] points at &TransMatrix_Cn[0][0]; the odd
        // part starts at row 1 of the 2-row table.
        let tr_even = &trans[0];
        let tr_odd = &trans[1];

        // Convert even elements of jac_coeff.
        for i in 0..=half_degree {
            let i_ptr_idx = i * max_m - (i + 1) * i / 2;
            let i_coeff_offset = double_dim * i;

            for idim in 1..=dimension {
                let mut value = 0.0f64;
                for j in i..=half_degree {
                    value += tr_even[i_ptr_idx + j] * jac_coeff[double_dim * j + idim - 1];
                }
                coefficients[i_coeff_offset + idim - 1] = value;
            }
        }

        if degree == 0 {
            return;
        }

        // Convert odd elements of jac_coeff.
        let half_degree_minus1 = (degree - 1) / 2;
        for i in 0..=half_degree_minus1 {
            let i_ptr_idx = i * max_m - (i + 1) * i / 2;
            let i_base_idx = (2 * i + 1) * dimension;
            let j_base_idx = (2 * i + 1) * dimension;

            for idim in 1..=dimension {
                let mut value = 0.0f64;
                let mut jj = j_base_idx + idim - 1;
                for j in i..=half_degree_minus1 {
                    value += tr_odd[i_ptr_idx + j] * jac_coeff[jj];
                    jj += double_dim;
                }
                coefficients[i_base_idx + idim - 1] = value;
            }
        }
    }
}
