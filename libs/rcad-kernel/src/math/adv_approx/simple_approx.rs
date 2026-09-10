//! OCCT AdvApprox_SimpleApprox (TKG3d).
//!
//! 1:1 translation of `AdvApprox_SimpleApprox.cxx` — approximates a function
//! on [First, Last] by a simple polynomial of lowest possible degree within
//! tolerance, using Gauss integration in a Jacobi basis.

use super::evaluator_function::EvaluatorFunction;
use crate::math::p_lib_jacobi::JacobiPolynomial;
use crate::math::plib::{eval_polynomial_flat, hermite_interpolate};
use crate::math::{GeomAbsShape, MatD};

/// OCCT AdvApprox_SimpleApprox.
pub struct SimpleApprox {
    total_num_ss: usize,
    total_dimension: usize,
    nb_gauss_points: usize,
    work_degree: usize,
    niv_constr: usize,
    jac_pol: JacobiPolynomial,
    tab_points: Vec<f64>,
    tab_weights: MatD,
    degree: usize,
    coeff: Vec<f64>,
    /// OCCT mySomTab (HArray1 with lower bound 0).
    som_tab: Vec<f64>,
    /// OCCT myDifTab (HArray1 with lower bound 0).
    dif_tab: Vec<f64>,
    max_error: Vec<f64>,
    average_error: Vec<f64>,
    done: bool,
}

impl SimpleApprox {
    /// OCCT AdvApprox_SimpleApprox ctor.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        total_dimension: usize,
        total_num_ss: usize,
        continuity: GeomAbsShape,
        work_degree: usize,
        nb_gauss_points: usize,
        jacobi_base: &JacobiPolynomial,
    ) -> Self {
        let niv_constr = match continuity {
            GeomAbsShape::C0 => 0,
            GeomAbsShape::C1 => 1,
            GeomAbsShape::C2 => 2,
            _ => panic!("Invalid Continuity"),
        };

        let degree_q = work_degree as i64 - 2 * (niv_constr as i64 + 1);

        // the extraction of the Legendre roots: TabPoints (0, NbGaussPoints/2).
        let mut tab_points = vec![0.0f64; nb_gauss_points / 2 + 1];
        jacobi_base.points(nb_gauss_points, &mut tab_points);

        // the extraction of the Gauss Weights: (0, NbGaussPoints/2) x
        // (0, DegreeQ).
        let rows = nb_gauss_points / 2 + 1;
        let cols = (degree_q + 1).max(0) as usize;
        let mut tab_weights = MatD::new(rows, cols);
        jacobi_base.weights(nb_gauss_points, &mut tab_weights);

        let coeff = vec![0.0f64; (work_degree + 1) * total_dimension];
        // OCCT mySomTab / myDifTab: HArray1(0, (NbGaussPoints/2 + 1) * dim - 1).
        let scratch = vec![0.0f64; (nb_gauss_points / 2 + 1) * total_dimension];

        SimpleApprox {
            total_num_ss,
            total_dimension,
            nb_gauss_points,
            work_degree,
            niv_constr,
            jac_pol: jacobi_base.clone(),
            tab_points,
            tab_weights,
            degree: 0,
            coeff,
            som_tab: scratch.clone(),
            dif_tab: scratch,
            max_error: Vec::new(),
            average_error: Vec::new(),
            done: false,
        }
    }

    /// OCCT AdvApprox_SimpleApprox::Perform.
    pub fn perform(
        &mut self,
        evaluator: &mut dyn EvaluatorFunction,
        local_dimension: &[i32],
        local_tolerances_array: &[f64],
        first: f64,
        last: f64,
        max_degree: usize,
    ) {
        self.done = false;
        let dimension = self.total_dimension;

        // ===== the computation of Rr(t) (the first part of Pp) =====
        let degree_r = 2 * self.niv_constr as i32 + 1;
        let _degree_q = self.work_degree as i32 - 2 * (self.niv_constr as i32 + 1);

        let first_last = [first, last];
        let mut result = vec![0.0f64; dimension];
        let fact = (last - first) / 2.0;

        // Constraints at First / Last (1..dimension x 0..nivConstr).
        let mut first_constr = MatD::new(dimension, self.niv_constr + 1);
        let mut last_constr = MatD::new(dimension, self.niv_constr + 1);

        for derive in (0..=self.niv_constr as i32).rev() {
            let err = evaluator.evaluate(&first_last, first, derive, &mut result);
            if err != 0 {
                return; // Evaluation error
            }
            if derive >= 1 {
                for v in result.iter_mut() {
                    *v *= fact;
                }
            }
            if derive == 2 {
                for v in result.iter_mut() {
                    *v *= fact;
                }
            }
            for idim in 1..=dimension {
                first_constr.set(idim, (derive + 1) as usize, result[idim - 1]);
            }
        }

        for derive in (0..=self.niv_constr as i32).rev() {
            let err = evaluator.evaluate(&first_last, last, derive, &mut result);
            if err != 0 {
                return; // Evaluation error
            }
            if derive >= 1 {
                for v in result.iter_mut() {
                    *v *= fact;
                }
            }
            if derive == 2 {
                for v in result.iter_mut() {
                    *v *= fact;
                }
            }
            for idim in 1..=dimension {
                last_constr.set(idim, (derive + 1) as usize, result[idim - 1]);
            }
        }

        hermite_interpolate(
            dimension,
            -1.0,
            1.0,
            self.niv_constr,
            self.niv_constr,
            &first_constr,
            &last_constr,
            &mut self.coeff,
        );

        // ===== the computation of the coefficients of Qq(t) =====
        let mut fti = vec![0.0f64; dimension];
        let mut rpti = vec![0.0f64; dimension];
        let mut rmti = vec![0.0f64; dimension];

        let alin = (last - first) / 2.0;
        let blin = (last + first) / 2.0;

        let dim_i = dimension as i32;
        // OCCT L176: i_idim = myTotalDimension — the paired-point slots start
        // at TotalDimension (slots 0..TotalDimension-1 are reserved for the
        // odd-NbGaussPoints case, L228-230).
        let mut i_idim = self.total_dimension as i32;
        for i in 1..=self.nb_gauss_points / 2 {
            let ti = self.tab_points[i];
            let tip = alin * ti + blin;
            let err = evaluator.evaluate(&first_last, tip, 0, &mut fti);
            if err != 0 {
                return; // Evaluation error
            }
            for idim in 1..=dim_i {
                // OCCT: mySomTab->SetValue(i_idim, ...) on an HArray1 with
                // lower bound 0 — the 0-based slot equals i_idim.
                let slot = i_idim as usize;
                self.som_tab[slot] = fti[idim as usize - 1];
                self.dif_tab[slot] = fti[idim as usize - 1];
                i_idim += 1;
            }
        }
        let mut i_idim = self.total_dimension as i32;
        for i in 1..=self.nb_gauss_points / 2 {
            let ti = self.tab_points[i];
            let tin = -alin * ti + blin;
            let err = evaluator.evaluate(&first_last, tin, 0, &mut fti);
            if err != 0 {
                return; // Evaluation error
            }
            eval_polynomial_flat(ti, 0, degree_r, dim_i, &self.coeff, &mut rpti);
            let ti_neg = -ti;
            eval_polynomial_flat(ti_neg, 0, degree_r, dim_i, &self.coeff, &mut rmti);

            for idim in 1..=dim_i {
                let slot = i_idim as usize;
                self.som_tab[slot] += fti[idim as usize - 1] - rpti[idim as usize - 1] - rmti[idim as usize - 1];
                self.dif_tab[slot] -= fti[idim as usize - 1] + rpti[idim as usize - 1] - rmti[idim as usize - 1];
                i_idim += 1;
            }
        }

        // for odd NbGaussPoints — the computation of [ F(0) - R(0) ].
        if self.nb_gauss_points % 2 == 1 {
            let ti = self.tab_points[0];
            let tip = blin;
            let err = evaluator.evaluate(&first_last, tip, 0, &mut fti);
            if err != 0 {
                return; // Evaluation error
            }
            eval_polynomial_flat(ti, 0, degree_r, dim_i, &self.coeff, &mut rpti);
            for idim in 1..=dim_i {
                self.som_tab[idim as usize - 1] = fti[idim as usize - 1] - rpti[idim as usize - 1];
                self.dif_tab[idim as usize - 1] = fti[idim as usize - 1] - rpti[idim as usize - 1];
            }
        }

        // the computation of Qq(t).  NOTE (OCCT verbatim): `Sum` is declared
        // once and NOT reset inside the odd-GaussPoints branch, so it carries
        // its previous value across iterations there.
        let mut sum = 0.0f64;
        let degree_q = self.work_degree as i32 - 2 * (self.niv_constr as i32 + 1);
        let mut k = 0i32;
        while k <= degree_q {
            for idim in 1..=dim_i {
                sum = 0.0;
                for i in 1..=self.nb_gauss_points / 2 {
                    sum += self.tab_weights.get(i + 1, (k + 1) as usize)
                        * self.som_tab[i * dimension + idim as usize - 1];
                }
                self.coeff[((k + degree_r + 1) * dim_i + idim - 1) as usize] = sum;
            }
            k += 2;
        }
        let mut k = 1i32;
        while k <= degree_q {
            for idim in 1..=dim_i {
                sum = 0.0;
                for i in 1..=self.nb_gauss_points / 2 {
                    sum += self.tab_weights.get(i + 1, (k + 1) as usize)
                        * self.dif_tab[i * dimension + idim as usize - 1];
                }
                self.coeff[((k + degree_r + 1) * dim_i + idim - 1) as usize] = sum;
            }
            k += 2;
        }
        if self.nb_gauss_points % 2 == 1 {
            for idim in 1..=dim_i {
                let mut k = 0i32;
                while k <= degree_q {
                    sum += self.tab_weights.get(1, (k + 1) as usize)
                        * self.som_tab[idim as usize - 1];
                    self.coeff[((k + degree_r + 1) * dim_i + idim - 1) as usize] = sum;
                    k += 2;
                }
            }
        }

        // the computing of NewDegree.
        let mut jac_coeff = vec![0.0f64; dimension * (self.work_degree + 1)];

        let mut max_err;
        let mut average_err;
        let mut new_degree_max = 0usize;
        let mut new_degree = 0usize;
        let mut max_error_out = 0.0f64;

        self.max_error = vec![0.0f64; self.total_num_ss];
        self.average_error = vec![0.0f64; self.total_num_ss];
        let mut rang_ss = 0usize;
        let mut rang_jac_coeff = 0usize;
        for numss in 1..=self.total_num_ss {
            let dim = local_dimension[numss - 1] as usize;
            let mut rang_coeff = 0usize;
            let mut rang_dim = 0usize;
            for k in 0..=self.work_degree {
                for idim in 1..=dim {
                    jac_coeff[rang_jac_coeff + rang_dim + idim - 1] =
                        self.coeff[rang_coeff + rang_ss + idim - 1];
                }
                rang_dim += dim;
                rang_coeff += dimension;
            }

            let jac_ss = &jac_coeff[rang_jac_coeff..];
            self.jac_pol.reduce_degree(
                dim,
                max_degree,
                local_tolerances_array[numss - 1],
                jac_ss,
                &mut new_degree,
                &mut max_error_out,
            );
            if new_degree > new_degree_max {
                new_degree_max = new_degree;
            }
            rang_ss += dim;
            rang_jac_coeff += (self.work_degree + 1) * dim;
        }

        // the computing of MaxError and AverageError.
        let mut rang_ss = 0usize;
        let mut rang_jac_coeff = 0usize;
        for numss in 1..=self.total_num_ss {
            let dim = local_dimension[numss - 1] as usize;
            let jac_ss = &jac_coeff[rang_jac_coeff..];
            max_err = self.jac_pol.max_error(dim, jac_ss, new_degree_max);
            self.max_error[numss - 1] = max_err;
            average_err = self.jac_pol.average_error(dim, jac_ss, new_degree_max);
            self.average_error[numss - 1] = average_err;
            rang_ss += dim;
            rang_jac_coeff += (self.work_degree + 1) * dim;
        }

        self.degree = new_degree_max;
        self.done = true;
    }

    /// OCCT AdvApprox_SimpleApprox::IsDone.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT AdvApprox_SimpleApprox::Degree.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// OCCT AdvApprox_SimpleApprox::Coefficients.
    pub fn coefficients(&self) -> &[f64] {
        &self.coeff
    }

    /// OCCT AdvApprox_SimpleApprox::MaxError(Index).
    pub fn max_error(&self, index: usize) -> f64 {
        self.max_error[index - 1]
    }

    /// OCCT AdvApprox_SimpleApprox::AverageError(Index).
    pub fn average_error(&self, index: usize) -> f64 {
        self.average_error[index - 1]
    }
}
