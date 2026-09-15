//! OCCT PLib_JacobiPolynomial (FoundationClasses/TKMath/PLib).
//!
//! 1:1 translation of `PLib_JacobiPolynomial.hxx` (L17-170) and
//! `PLib_JacobiPolynomial.cxx` (L1-553, including the D0123 basis recursion,
//! FORTRAN MPOJAC port, which the older partial `math::p_lib_jacobi` module
//! does not carry).
//!
//! Container mapping: OCCT `NCollection_Array1<double>` basis/coefficient
//! arrays are passed as `&mut [f64]` slices; every OCCT call site of these
//! APIs uses 0-based arrays (Lower() == 0), so slice indexing equals the
//! OCCT indexing.  2D outputs (`Weights`) and internal `math_Matrix` /
//! `math_Vector` locals use the bound-aware `math_matrix::{Matrix, Vector}`
//! mirrors.  The `double& theJacCoeff` coefficient arrays are flat
//! `[f64]` slices addressed exactly like the OCCT pointer arithmetic.

use super::super::gauss_points::gauss_points;
use super::super::math_matrix::{Matrix, Vector};
use super::super::p_lib_jacobi_data::{
    MAX_VALUES_DB_C0, MAX_VALUES_DB_C1, MAX_VALUES_DB_C2, TRANS_MATRIX_C0, TRANS_MATRIX_C1,
    TRANS_MATRIX_C2, WEIGHTS_DB0_C0, WEIGHTS_DB0_C1, WEIGHTS_DB0_C2, WEIGHTS_DB_C0, WEIGHTS_DB_C1,
    WEIGHTS_DB_C2,
};
use super::super::GeomAbsShape;
use super::super::p_lib_jacobi::niv_constr;
use super::jacobi_coefficients_data::get_jacobi_coefficients;

/// The possible values for NbGaussPoints (PLib_JacobiPolynomial.cxx L32-34).
const THE_NB_GAUSS_POINTS_8: i32 = 8;
const THE_NB_GAUSS_POINTS_10: i32 = 10;
const THE_NB_GAUSS_POINTS_15: i32 = 15;
const THE_NB_GAUSS_POINTS_20: i32 = 20;
const THE_NB_GAUSS_POINTS_25: i32 = 25;
const THE_NB_GAUSS_POINTS_30: i32 = 30;
const THE_NB_GAUSS_POINTS_40: i32 = 40;
const THE_NB_GAUSS_POINTS_50: i32 = 50;
const THE_NB_GAUSS_POINTS_61: i32 = 61;

/// OCCT THE_INVALID_VALUE (PLib_JacobiPolynomial.cxx L36).
const THE_INVALID_VALUE: f64 = -999.0;

/// Maximum supported polynomial degree (PLib_JacobiPolynomial.cxx L39).
const THE_MAX_DEGREE: i32 = 30;

/// Lookup tables for database pointers indexed by myNivConstr
/// (PLib_JacobiPolynomial.cxx L42-47).
const THE_WEIGHTS_DB: [&[f64]; 3] = [&WEIGHTS_DB_C0, &WEIGHTS_DB_C1, &WEIGHTS_DB_C2];
const THE_WEIGHTS_DB0: [&[f64]; 3] = [&WEIGHTS_DB0_C0, &WEIGHTS_DB0_C1, &WEIGHTS_DB0_C2];
const THE_MAX_VALUES_DB: [&[f64]; 3] = [&MAX_VALUES_DB_C0, &MAX_VALUES_DB_C1, &MAX_VALUES_DB_C2];

/// OCCT PLib_JacobiPolynomial - Jacobi polynomial basis relative to an order
/// of constraint (PLib_JacobiPolynomial.hxx L27-52):
/// P(t) = R(t) + W(t) * Q(t), W(t) = (1-t^2)^(2*nivConstr+2).
#[derive(Debug, Clone)]
pub struct JacobiPolynomial {
    work_degree: i32,
    niv_constr: i32,
    degree: i32,
}

impl JacobiPolynomial {
    /// OCCT PLib_JacobiPolynomial::PLib_JacobiPolynomial
    /// (PLib_JacobiPolynomial.cxx L52-66).
    pub fn new(the_work_degree: i32, the_constraint_order: GeomAbsShape) -> Self {
        // OCCT PLib::NivConstr returns int; the rcad helper returns usize.
        let my_niv_constr = niv_constr(the_constraint_order) as i32;
        let my_degree = the_work_degree - 2 * (my_niv_constr + 1);
        assert!(
            my_degree >= 0,
            "Standard_ConstructionError: WorkDegree too small for given ConstraintOrder"
        );
        assert!(
            my_degree <= THE_MAX_DEGREE,
            "Standard_ConstructionError: Invalid Degree"
        );
        JacobiPolynomial {
            work_degree: the_work_degree,
            niv_constr: my_niv_constr,
            degree: my_degree,
        }
    }

    /// OCCT PLib_JacobiPolynomial::WorkDegree (hxx L149).
    #[inline]
    pub fn work_degree(&self) -> i32 {
        self.work_degree
    }

    /// OCCT PLib_JacobiPolynomial::NivConstr (hxx L152).
    #[inline]
    pub fn niv_constr(&self) -> i32 {
        self.niv_constr
    }

    /// OCCT PLib_JacobiPolynomial::Points (PLib_JacobiPolynomial.cxx L70-100).
    /// `the_tab_points` is the caller-allocated array over (0, NbGaussPoints/2).
    pub fn points(&self, the_nb_gauss_points: i32, the_tab_points: &mut Vector) {
        if (the_nb_gauss_points != THE_NB_GAUSS_POINTS_8
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_10
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_15
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_20
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_25
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_30
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_40
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_50
            && the_nb_gauss_points != THE_NB_GAUSS_POINTS_61)
            || the_nb_gauss_points <= self.degree
        {
            panic!("Standard_ConstructionError: Invalid NbGaussPoints");
        }

        let mut a_decreasing_points = Vector::new(1, the_nb_gauss_points);

        gauss_points(
            the_nb_gauss_points as usize,
            &mut a_decreasing_points.data,
        );

        // theTabPoints consist of only positive increasing values
        for i in 1..=the_nb_gauss_points / 2 {
            let v = a_decreasing_points.get(the_nb_gauss_points / 2 - i + 1);
            the_tab_points.set(i, v);
        }
        if the_nb_gauss_points % 2 == 1 {
            the_tab_points.set(0, 0.0);
        } else {
            the_tab_points.set(0, THE_INVALID_VALUE);
        }
    }

    /// OCCT PLib_JacobiPolynomial::Weights (PLib_JacobiPolynomial.cxx
    /// L104-186).  `the_tab_weights` is the caller-allocated matrix over
    /// (0, NbGaussPoints/2) x (0, Degree).
    pub fn weights(&self, the_nb_gauss_points: i32, the_tab_weights: &mut Matrix) {
        // Zero-initialize entire output array to ensure columns beyond myDegree
        // (which are not populated below) contain valid finite values.
        the_tab_weights.init(0.0);

        let mut a_db_pointer: usize = 0;
        let a_min_degree = 2 * (self.niv_constr + 1);

        // Calculate offset into the weights database
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_8 {
            a_db_pointer += (THE_NB_GAUSS_POINTS_8 * (THE_NB_GAUSS_POINTS_8 - a_min_degree) / 2) as usize;
        }
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_10 {
            a_db_pointer += (THE_NB_GAUSS_POINTS_10 * (THE_NB_GAUSS_POINTS_10 - a_min_degree) / 2) as usize;
        }
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_15 {
            a_db_pointer += (((THE_NB_GAUSS_POINTS_15 - 1) / 2) * (THE_NB_GAUSS_POINTS_15 - a_min_degree)) as usize;
        }
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_20 {
            a_db_pointer += (THE_NB_GAUSS_POINTS_20 * (THE_NB_GAUSS_POINTS_20 - a_min_degree) / 2) as usize;
        }
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_25 {
            a_db_pointer += (((THE_NB_GAUSS_POINTS_25 - 1) / 2) * (THE_NB_GAUSS_POINTS_25 - a_min_degree)) as usize;
        }
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_30 {
            a_db_pointer += (THE_NB_GAUSS_POINTS_30 * (THE_NB_GAUSS_POINTS_30 - a_min_degree) / 2) as usize;
        }
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_40 {
            a_db_pointer += (THE_NB_GAUSS_POINTS_40 * (THE_NB_GAUSS_POINTS_40 - a_min_degree) / 2) as usize;
        }
        if the_nb_gauss_points > THE_NB_GAUSS_POINTS_50 {
            a_db_pointer += (THE_NB_GAUSS_POINTS_50 * (THE_NB_GAUSS_POINTS_50 - a_min_degree) / 2) as usize;
        }
        let a_db = &THE_WEIGHTS_DB[self.niv_constr as usize][a_db_pointer..];

        // Copy TabWeightsDB into theTabWeights (explicit loops for safe 2D array access)
        let mut a_db_index = 0usize;
        let a_half_points = the_nb_gauss_points / 2;
        for j in 0..=self.degree {
            for i in 1..=a_half_points {
                the_tab_weights.set(i, j, a_db[a_db_index]);
                a_db_index += 1;
            }
        }

        if the_nb_gauss_points % 2 == 1 {
            // theNbGaussPoints is odd - fill row 0 with special values
            let mut a_db_pointer0: usize = 0;
            let a_db0 = THE_WEIGHTS_DB0[self.niv_constr as usize];

            if the_nb_gauss_points > THE_NB_GAUSS_POINTS_15 {
                a_db_pointer0 += ((THE_NB_GAUSS_POINTS_15 - 1 - a_min_degree) / 2 + 1) as usize;
            }
            if the_nb_gauss_points > THE_NB_GAUSS_POINTS_25 {
                a_db_pointer0 += ((THE_NB_GAUSS_POINTS_25 - 1 - a_min_degree) / 2 + 1) as usize;
            }

            // Overwrite even columns of row 0 with data from database
            let mut j = 0;
            while j <= self.degree {
                the_tab_weights.set(0, j, a_db0[a_db_pointer0]);
                a_db_pointer0 += 1;
                j += 2;
            }
        } else {
            // Fill row 0 with THE_INVALID_VALUE for even theNbGaussPoints
            // (explicit column loop)
            for j in 0..=self.degree {
                the_tab_weights.set(0, j, THE_INVALID_VALUE);
            }
        }
    }

    /// OCCT PLib_JacobiPolynomial::MaxValue (PLib_JacobiPolynomial.cxx
    /// L190-197).  `the_tab_max` is the caller-allocated array over
    /// (0, myDegree + 1).
    pub fn max_value(&self, the_tab_max: &mut [f64]) {
        let a_db_pointer = THE_MAX_VALUES_DB[self.niv_constr as usize];
        for (i, v) in the_tab_max.iter_mut().enumerate() {
            *v = a_db_pointer[i];
        }
    }

    /// OCCT PLib_JacobiPolynomial::MaxError (PLib_JacobiPolynomial.cxx
    /// L201-227).
    pub fn max_error(&self, the_dimension: i32, the_jac_coeff: &[f64], the_new_degree: i32) -> f64 {
        // Buffering on stack to avoid dynamic allocation in this frequently
        // called method
        let mut a_tab_max_buf = vec![0.0; (self.degree + 2) as usize];
        self.max_value(&mut a_tab_max_buf);
        let a_tab_max = &a_tab_max_buf[..];

        let a_beg_idx = 2 * (self.niv_constr + 1);
        let a_cut_idx = a_beg_idx.max(the_new_degree + 1);

        let mut a_max_err_dim = Vector::new_init(1, the_dimension, 0.0);

        for a_dim_idx in 1..=the_dimension {
            for a_coeff_idx in a_cut_idx..=self.work_degree {
                let a_coeff_value = the_jac_coeff[(a_coeff_idx * the_dimension + a_dim_idx - 1) as usize];
                let a_basis_max = a_tab_max[(a_coeff_idx - a_beg_idx) as usize];
                let v = a_max_err_dim.get(a_dim_idx) + a_coeff_value.abs() * a_basis_max;
                a_max_err_dim.set(a_dim_idx, v);
            }
        }

        // aMaxErrDim.Norm()
        a_max_err_dim.data.v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    /// OCCT PLib_JacobiPolynomial::ReduceDegree (PLib_JacobiPolynomial.cxx
    /// L231-289).
    pub fn reduce_degree(
        &self,
        the_dimension: i32,
        the_max_degree: i32,
        the_tol: f64,
        the_jac_coeff: &[f64],
        the_new_degree: &mut i32,
        the_max_error: &mut f64,
    ) {
        let an_idx = 2 * (self.niv_constr + 1) - 1;
        let a_cut_idx = an_idx + 1;

        let mut a_max_err_dim = Vector::new_init(1, the_dimension, 0.0);
        let mut a_tab_max = vec![0.0; (self.degree + 2) as usize];
        self.max_value(&mut a_tab_max);

        *the_new_degree = an_idx;
        *the_max_error = 0.0;

        // Search for theNewDegree from high degree to low
        for i in (a_cut_idx..=self.work_degree).rev() {
            // Accumulate error contribution for all dimensions
            let i_offset = (i * the_dimension) as usize;
            for idim in 1..=the_dimension {
                let v = a_max_err_dim.get(idim)
                    + the_jac_coeff[i_offset + (idim - 1) as usize].abs()
                        * a_tab_max[(i - a_cut_idx) as usize];
                a_max_err_dim.set(idim, v);
            }

            // aMaxErrDim.Norm()
            let an_error = a_max_err_dim.data.v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if an_error > the_tol && i <= the_max_degree {
                *the_new_degree = i;
                break;
            }
            *the_max_error = an_error;
        }

        // Fallback: find last non-negligible coefficient
        if *the_new_degree == an_idx {
            const AN_EPS: f64 = 1.0e-9;
            *the_new_degree = 0;

            for i in (1..=an_idx).rev() {
                let mut a_bid = 0.0;
                let i_offset = (i * the_dimension) as usize;
                for idim in 1..=the_dimension {
                    a_bid += the_jac_coeff[i_offset + (idim - 1) as usize].abs();
                }
                if a_bid > AN_EPS {
                    *the_new_degree = i;
                    break;
                }
            }
        }
    }

    /// OCCT PLib_JacobiPolynomial::AverageError (PLib_JacobiPolynomial.cxx
    /// L293-312).
    pub fn average_error(&self, the_dimension: i32, the_jac_coeff: &[f64], the_new_degree: i32) -> f64 {
        let a_cut_idx = (2 * (self.niv_constr + 1) + 1).max(the_new_degree + 1);
        let mut an_average_err = 0.0;

        // Compute sum of squares of coefficients beyond theNewDegree
        for idim in 1..=the_dimension {
            for i in a_cut_idx..=self.degree {
                let a_jac_coeff = the_jac_coeff[(i * the_dimension + idim - 1) as usize];
                an_average_err += a_jac_coeff * a_jac_coeff;
            }
        }

        (an_average_err / 2.0).sqrt()
    }

    /// OCCT PLib_JacobiPolynomial::ToCoefficients
    /// (PLib_JacobiPolynomial.cxx L316-372).  Converts the polynomial
    /// P(t) = R(t) + W(t) Q(t) in the canonical base.  Both coefficient
    /// arrays use the Lower() == 0 slice convention.
    pub fn to_coefficients(
        &self,
        the_dimension: i32,
        the_degree: i32,
        the_jac_coeff: &[f64],
        the_coefficients: &mut [f64],
    ) {
        const A_MAX_M: i32 = THE_MAX_DEGREE + 1;
        let a_half_degree = the_degree / 2;
        let a_double_dim = 2 * the_dimension;

        // OCCT THE_TRANS_MATRIX[niv] is one flat block of 2 * (A_MAX_M
        // choose 2) doubles; the rcad TRANS_MATRIX_Cx tables store it as the
        // two halves [f64; A_MAX_M*(A_MAX_M+1)/2], so the odd-phase pointer
        // bump `aTrPointer += A_MAX_M*(A_MAX_M+1)/2` selects the second half.
        let trans_matrix: [&[f64]; 2] = match self.niv_constr {
            0 => [&TRANS_MATRIX_C0[0], &TRANS_MATRIX_C0[1]],
            1 => [&TRANS_MATRIX_C1[0], &TRANS_MATRIX_C1[1]],
            _ => [&TRANS_MATRIX_C2[0], &TRANS_MATRIX_C2[1]],
        };
        let mut a_tr_phase = 0usize;
        let mut a_tr_pointer: &[f64] = trans_matrix[a_tr_phase];

        // Convert even elements of theJacCoeff
        for i in 0..=a_half_degree {
            let i_ptr_idx = (i * A_MAX_M - (i + 1) * i / 2) as usize;
            let i_coeff_offset = (a_double_dim * i) as usize;

            for idim in 1..=the_dimension {
                let mut a_value = 0.0;
                for j in i..=a_half_degree {
                    a_value += a_tr_pointer[i_ptr_idx + j as usize]
                        * the_jac_coeff[(a_double_dim * j + idim - 1) as usize];
                }
                the_coefficients[i_coeff_offset + (idim - 1) as usize] = a_value;
            }
        }

        if the_degree == 0 {
            return;
        }

        // Convert odd elements of theJacCoeff
        a_tr_phase = 1;
        a_tr_pointer = trans_matrix[a_tr_phase];
        let a_half_degree_minus1 = (the_degree - 1) / 2;

        for i in 0..=a_half_degree_minus1 {
            let i_ptr_idx = (i * A_MAX_M - (i + 1) * i / 2) as usize;
            let i_base_idx = ((2 * i + 1) * the_dimension) as usize;
            let j_base_idx = ((2 * i + 1) * the_dimension) as usize;

            for idim in 1..=the_dimension {
                let mut a_value = 0.0;
                let mut jj = j_base_idx + (idim - 1) as usize;

                for j in i..=a_half_degree_minus1 {
                    a_value += a_tr_pointer[i_ptr_idx + j as usize] * the_jac_coeff[jj];
                    jj += a_double_dim as usize;
                }
                the_coefficients[i_base_idx + (idim - 1) as usize] = a_value;
            }
        }
    }

    /// OCCT PLib_JacobiPolynomial::D0123 (PLib_JacobiPolynomial.cxx L379-515),
    /// common part of D0, D1, D2, D3 (FORTRAN MPOJAC port).
    ///
    /// Architecture deviation: OCCT passes the SAME array several times
    /// (`D0: D0123(0, U, V, V, V, V)`), which the C++ aliasing rules allow
    /// and Rust `&mut` forbids.  The derivative outputs therefore arrive as
    /// `Option`, consumed only for the derivative order actually requested;
    /// the public D0..D3 wrappers below reproduce the OCCT call forms.
    fn d0123(
        &self,
        the_n_deriv: i32,
        the_u: f64,
        the_basis_value: &mut [f64],
        the_basis_d1: Option<&mut [f64]>,
        the_basis_d2: Option<&mut [f64]>,
        the_basis_d3: Option<&mut [f64]>,
    ) {
        let a_hermit_niv_constr = 2 * (self.niv_constr + 1);

        // Consume the derivative outputs for the requested order (see the
        // Option-based signature note above); the dummy fallbacks are never
        // dereferenced because every use below is guarded by theNDeriv.
        let mut empty1: Vec<f64> = Vec::new();
        let mut empty2: Vec<f64> = Vec::new();
        let mut empty3: Vec<f64> = Vec::new();
        let mut the_basis_d1 = the_basis_d1;
        let mut the_basis_d2 = the_basis_d2;
        let mut the_basis_d3 = the_basis_d3;
        let the_basis_d1: &mut [f64] = if the_n_deriv >= 1 {
            the_basis_d1.take().expect("BasisD1 required for NDeriv >= 1")
        } else {
            &mut empty1[..]
        };
        let the_basis_d2: &mut [f64] = if the_n_deriv >= 2 {
            the_basis_d2.take().expect("BasisD2 required for NDeriv >= 2")
        } else {
            &mut empty2[..]
        };
        let the_basis_d3: &mut [f64] = if the_n_deriv >= 3 {
            the_basis_d3.take().expect("BasisD3 required for NDeriv == 3")
        } else {
            &mut empty3[..]
        };

        // Get pre-computed coefficients from static cache (zero per-instance
        // overhead!)
        let (a_t_norm, a_cof_a, a_cof_b, a_denom) =
            get_jacobi_coefficients(self.niv_constr as usize, self.degree as usize);

        // --- Positionements triviaux -----
        // OCCT reads BasisValue.Lower() etc.; all call sites pass 0-based
        // arrays, so aBeg0 = aBeg1 = aBeg2 = aBeg3 = 0.
        let a_beg0 = 0usize;
        let a_beg1 = 0usize;
        let a_beg2 = 0usize;
        let a_beg3 = 0usize;

        if self.degree == 0 {
            the_basis_value[a_beg0] = 1.0;
            if the_n_deriv >= 1 {
                the_basis_d1[a_beg1] = 0.0;
                if the_n_deriv >= 2 {
                    the_basis_d2[a_beg2] = 0.0;
                    if the_n_deriv == 3 {
                        the_basis_d3[a_beg3] = 0.0;
                    }
                }
            }
        } else {
            the_basis_value[a_beg0] = 1.0;
            let an_aux = (a_hermit_niv_constr + 1) as f64;
            the_basis_value[a_beg0 + 1] = an_aux * the_u;
            if the_n_deriv >= 1 {
                the_basis_d1[a_beg1] = 0.0;
                the_basis_d1[a_beg1 + 1] = an_aux;
                if the_n_deriv >= 2 {
                    the_basis_d2[a_beg2] = 0.0;
                    the_basis_d2[a_beg2 + 1] = 0.0;
                    if the_n_deriv == 3 {
                        the_basis_d3[a_beg3] = 0.0;
                        the_basis_d3[a_beg3 + 1] = 0.0;
                    }
                }
            }
        }

        // --- Positionement par reccurence
        if self.degree > 1 {
            if the_n_deriv == 0 {
                for i in 2..=self.degree {
                    let iu = i as usize;
                    the_basis_value[a_beg0 + iu] = (a_cof_a[iu] * the_u
                        * the_basis_value[a_beg0 + iu - 1]
                        + a_cof_b[iu] * the_basis_value[a_beg0 + iu - 2])
                        * a_denom[iu];
                }
            } else {
                for i in 2..=self.degree {
                    let i0 = (i + a_beg0 as i32) as usize;
                    let i1 = (i + a_beg1 as i32) as usize;
                    let a_cof_a = a_cof_a[i as usize];
                    let a_cof_b = a_cof_b[i as usize];
                    let a_denom = a_denom[i as usize];

                    the_basis_value[i0] = (a_cof_a * the_u * the_basis_value[i0 - 1]
                        + a_cof_b * the_basis_value[i0 - 2])
                        * a_denom;
                    the_basis_d1[i1] = (a_cof_a
                        * (the_u * the_basis_d1[i1 - 1] + the_basis_value[i0 - 1])
                        + a_cof_b * the_basis_d1[i1 - 2])
                        * a_denom;
                    if the_n_deriv >= 2 {
                        let i2 = (i + a_beg2 as i32) as usize;
                        the_basis_d2[i2] = (a_cof_a
                            * (the_u * the_basis_d2[i2 - 1] + 2.0 * the_basis_d1[i1 - 1])
                            + a_cof_b * the_basis_d2[i2 - 2])
                            * a_denom;
                        if the_n_deriv == 3 {
                            let i3 = (i + a_beg3 as i32) as usize;
                            the_basis_d3[i3] = (a_cof_a
                                * (the_u * the_basis_d3[i3 - 1] + 3.0 * the_basis_d2[i2 - 1])
                                + a_cof_b * the_basis_d3[i3 - 2])
                                * a_denom;
                        }
                    }
                }
            }
        }

        // Normalization
        if the_n_deriv == 0 {
            for i in 0..=self.degree {
                let iu = i as usize;
                the_basis_value[a_beg0 + iu] *= a_t_norm[iu];
            }
        } else {
            for i in 0..=self.degree {
                let iu = i as usize;
                let a_norm_value = a_t_norm[iu];
                the_basis_value[iu + a_beg0] *= a_norm_value;
                the_basis_d1[iu + a_beg1] *= a_norm_value;
                if the_n_deriv >= 2 {
                    the_basis_d2[iu + a_beg2] *= a_norm_value;
                    if the_n_deriv >= 3 {
                        the_basis_d3[iu + a_beg3] *= a_norm_value;
                    }
                }
            }
        }
    }

    /// OCCT PLib_JacobiPolynomial::D0 (PLib_JacobiPolynomial.cxx L519-522).
    pub fn d0(&self, the_u: f64, the_basis_value: &mut [f64]) {
        self.d0123(0, the_u, the_basis_value, None, None, None);
    }

    /// OCCT PLib_JacobiPolynomial::D1 (PLib_JacobiPolynomial.cxx L526-531).
    pub fn d1(&self, the_u: f64, the_basis_value: &mut [f64], the_basis_d1: &mut [f64]) {
        self.d0123(
            1,
            the_u,
            the_basis_value,
            Some(the_basis_d1),
            None,
            None,
        );
    }

    /// OCCT PLib_JacobiPolynomial::D2 (PLib_JacobiPolynomial.cxx L535-541).
    pub fn d2(
        &self,
        the_u: f64,
        the_basis_value: &mut [f64],
        the_basis_d1: &mut [f64],
        the_basis_d2: &mut [f64],
    ) {
        self.d0123(
            2,
            the_u,
            the_basis_value,
            Some(the_basis_d1),
            Some(the_basis_d2),
            None,
        );
    }

    /// OCCT PLib_JacobiPolynomial::D3 (PLib_JacobiPolynomial.cxx L545-552).
    pub fn d3(
        &self,
        the_u: f64,
        the_basis_value: &mut [f64],
        the_basis_d1: &mut [f64],
        the_basis_d2: &mut [f64],
        the_basis_d3: &mut [f64],
    ) {
        self.d0123(
            3,
            the_u,
            the_basis_value,
            Some(the_basis_d1),
            Some(the_basis_d2),
            Some(the_basis_d3),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT PLib_JacobiPolynomial D0 basis, hand-derived from the MPOJAC
    /// recurrence (cxx L379-515) for NivConstr=0, WorkDegree=4 => Degree=2:
    ///   J0 = 1,  J1 = 3u,  J2 = (CofA2*u*J1 + CofB2*J0)*Denom2
    ///      = (336u*3u - 144)*(1/144) = 7u^2 - 1,
    /// then normalization by TNorm.  The C0 table squares are the exact
    /// rationals TNorm^2 = [15/16, 35/48, 45/64] (hand-spotted), giving
    /// at u = 1/2:
    ///   J0 = sqrt(15/16) = 0.9682458365518543
    ///   J1 = (3/2)*sqrt(35/48) = 1.2808688457449498
    ///   J2 = (7/4 - 1)*sqrt(45/64) = 0.6288941186718159
    /// and D1 = [0, 3*sqrt(35/48), 14u*sqrt(45/64)] = [0, 2.5617376914898996,
    /// 5.869678440936948].
    #[test]
    fn jacobi_d0_d1_degree2_c0_hand_recurrence() {
        let jacobi = JacobiPolynomial::new(4, GeomAbsShape::C0);

        let u = 0.5;
        let mut v = [0.0f64; 3];
        let mut d1 = [0.0f64; 3];
        jacobi.d1(u, &mut v, &mut d1);

        let tn0 = (15.0f64 / 16.0).sqrt();
        let tn1 = (35.0f64 / 48.0).sqrt();
        let tn2 = (45.0f64 / 64.0).sqrt();

        assert!((v[0] - tn0).abs() < 1e-14, "v0 = {}", v[0]);
        assert!((v[1] - 1.5 * tn1).abs() < 1e-14, "v1 = {}", v[1]);
        assert!((v[2] - 0.75 * tn2).abs() < 1e-14, "v2 = {}", v[2]);

        assert!(d1[0].abs() < 1e-14);
        assert!((d1[1] - 3.0 * tn1).abs() < 1e-14, "d1[1] = {}", d1[1]);
        assert!((d1[2] - 7.0 * tn2).abs() < 1e-14, "d1[2] = {}", d1[2]);
    }

    /// OCCT PLib_JacobiPolynomial::ToCoefficients hand-check for Degree=2,
    /// C0 (Dimension = 1), traced over OCCT PLib_JacobiPolynomial.cxx
    /// L316-372 with the TransMatrix_C0 data of PLib_JacobiPolynomial_Data.pxx.
    /// The even phase reads row 0 entries T[0] = 1 (canonical[0] <- JacCoeff[0])
    /// and T[1] = +TN0 (canonical[0] <- JacCoeff[2]), then row 1 entries
    /// T[31] = -TN0 (canonical[2] <- JacCoeff[2]); the odd phase reads
    /// TransMatrix_C0[1][0] = 1 (canonical[1] <- JacCoeff[1]), with
    /// TN0 = 0.9682458365518542... = sqrt(15/16) (hand-spotted rational).
    /// Hence for JacCoeff = [1, 2, 3]:
    ///     canonical = [1 + 3*TN0, 2, -3*TN0].
    #[test]
    fn jacobi_to_coefficients_degree2_hand() {
        let jacobi = JacobiPolynomial::new(4, GeomAbsShape::C0);
        let tn0 = (15.0f64 / 16.0).sqrt();

        let jac_coeff = [1.0, 2.0, 3.0];
        let mut coefficients = [0.0f64; 3];
        jacobi.to_coefficients(1, 2, &jac_coeff, &mut coefficients);

        assert!((coefficients[0] - (1.0 + 3.0 * tn0)).abs() < 1e-14);
        assert!((coefficients[1] - 2.0).abs() < 1e-14);
        assert!((coefficients[2] - (-3.0 * tn0)).abs() < 1e-14);
    }
}
