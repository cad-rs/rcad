//! OCCT AdvApprox_ApproxAFunction (TKG3d).
//!
//! 1:1 translation of `AdvApprox_ApproxAFunction.cxx` (static
//! `PrepareConvert` L123-313, static `Approximation` L364-599, `Perform`
//! L663-956 and the accessors).

use super::cutting::Cutting;
use super::evaluator_function::EvaluatorFunction;
use super::simple_approx::SimpleApprox;
use super::continuity_order;
use crate::math::bspl_lib::nb_poles;
use crate::math::convert_comp_polynomial_to_poles::ConvertCompPolynomialToPoles;
use crate::math::p_lib_jacobi::JacobiPolynomial;
use crate::math::plib::jacobi_parameters;
use crate::math::GeomAbsShape;

/// Determine local continuities (OCCT static PrepareConvert, L123-313).
/// All arrays are flat with the OCCT index layouts preserved:
/// - `num_coeff_per_curve[1..=num_curves]`
/// - `coefficients[(icurve-1)*dimension*real_degree ... ]`
/// - `polynomial_intervals[(icurve, 1..2)]` row-major with 2 columns
/// - `true_intervals[1..=num_curves+1]`
/// - `error_max[1 .. num_curves*nbspace]`
#[allow(clippy::too_many_arguments)]
pub(crate) fn prepare_convert(
    num_curves: i32,
    max_degree: i32,
    continuity_order: i32,
    num1dss: i32,
    num2dss: i32,
    num3dss: i32,
    num_coeff_per_curve: &[i32],
    coefficients: &mut [f64],
    polynomial_intervals: &[f64],
    true_intervals: &[f64],
    local_tolerance: &[f64],
    error_max: &mut [f64],
    continuity: &mut [i32],
) {
    let dimension = (num1dss + 2 * num2dss + 3 * num3dss) as usize;
    let nbspace = (num1dss + num2dss + num3dss) as usize;
    let real_degree = (max_degree + 1).max(2 * continuity_order + 2) as usize;

    // Init.
    for v in continuity.iter_mut() {
        *v = 0;
    }
    if continuity_order == 0 {
        return;
    }

    for icurve in 1..num_curves {
        // Init and positioning at the node.
        let mut is_ci = true;
        let coef1 = (icurve as usize - 1) * dimension * real_degree;
        let coef2 = coef1 + dimension * real_degree;
        let deg1 = num_coeff_per_curve[icurve as usize - 1] - 1;
        let deg2 = num_coeff_per_curve[icurve as usize] - 1;

        // Result holds 2 * (ContinuityOrder+1) * Dimension values: Res1 =
        // first half, Res2 = second half.
        let half = (continuity_order as usize + 1) * dimension;
        let mut result = vec![0.0f64; 2 * half];
        {
            let (res1, res2) = result.split_at_mut(half);
            crate::math::plib::eval_polynomial_flat(
                polynomial_intervals[((icurve - 1) * 2 + 1) as usize],
                continuity_order,
                deg1,
                dimension as i32,
                &coefficients[coef1..],
                res1,
            );
            crate::math::plib::eval_polynomial_flat(
                polynomial_intervals[icurve as usize * 2],
                continuity_order,
                deg2,
                dimension as i32,
                &coefficients[coef2..],
                res2,
            );
        }

        // Check in each subspace.
        for iordre in 1..=continuity_order {
            if !is_ci {
                break;
            }
            // fixing a bug PRO18577.
            let toler = 1.0e-5;

            let f1_dividend =
                polynomial_intervals[((icurve - 1) * 2 + 1) as usize] - polynomial_intervals[((icurve - 1) * 2) as usize];
            let f2_dividend = polynomial_intervals[(icurve * 2 + 1) as usize] - polynomial_intervals[(icurve * 2) as usize];
            let f1_divizor = true_intervals[icurve as usize] - true_intervals[(icurve - 1) as usize];
            let f2_divizor = true_intervals[(icurve + 1 - 1) as usize + 1] - true_intervals[icurve as usize];
            let facteur1;
            let facteur2;
            if f1_divizor.abs() < toler {
                facteur1 = 0.0;
            } else {
                let fract1 = f1_dividend / f1_divizor;
                facteur1 = fract1.powi(iordre);
            }
            if f2_divizor.abs() < toler {
                facteur2 = 0.0;
            } else {
                let fract2 = f2_dividend / f2_divizor;
                facteur2 = fract2.powi(iordre);
            }
            let normal1 = f1_divizor.powi(iordre);
            let normal2 = f2_divizor.powi(iordre);

            // Prec / Suivant scratch (1..NbSpace).
            let mut prec = vec![0.0f64; nbspace + 1];
            let mut suivant = vec![0.0f64; nbspace + 1];

            let mut idim = 1usize;
            // Val1 = Res1 + iordre*Dimension (OCCT pointer arithmetics with
            // the -1 offset baked in: Val1[ii] == Res1[iordre*dim + ii - 1]).
            let base1 = iordre as usize * dimension;
            let base2 = iordre as usize * dimension;

            // 1D subspaces.
            for ii in 1..=num1dss {
                if !is_ci {
                    break;
                }
                let idx = idim; // LocalTolerance(idim), 1-based
                let eps = local_tolerance[idx - 1] * 0.01;
                let v1 = result[base1 + ii as usize - 1] * facteur1;
                let v2 = result[half + base2 + ii as usize - 1] * facteur2;
                let diff = (v1 - v2).abs();
                let moy = (v1 + v2).abs();
                // A first check on the relative value.
                if diff > moy * 1.0e-9 {
                    prec[idim] = diff * normal1;
                    suivant[idim] = diff * normal2;
                    // And a second check on the upper bound of the error.
                    if prec[idim] > eps || suivant[idim] > eps {
                        is_ci = false;
                    }
                } else {
                    prec[idim] = 0.0;
                    suivant[idim] = 0.0;
                }
                idim += 1;
            }

            // 2D subspaces (pairs of coordinates).
            for ii in 1..=num2dss {
                if !is_ci {
                    break;
                }
                let idx = idim;
                let eps = local_tolerance[idx - 1] * 0.01;
                let o1 = base1 + (num1dss as usize) + (ii as usize - 1) * 2;
                let o2 = half + base2 + (num1dss as usize) + (ii as usize - 1) * 2;
                let v1x = result[o1] * facteur1;
                let v1y = result[o1 + 1] * facteur1;
                let v2x = result[o2] * facteur2;
                let v2y = result[o2 + 1] * facteur2;
                let diff = (v1x - v2x).abs() + (v1y - v2y).abs();
                let moy = (v1x + v2x).abs() + (v1y + v2y).abs();
                if diff > moy * 1.0e-9 {
                    prec[idim] = diff * normal1;
                    suivant[idim] = diff * normal2;
                    if prec[idim] > eps || suivant[idim] > eps {
                        is_ci = false;
                    }
                } else {
                    prec[idim] = 0.0;
                    suivant[idim] = 0.0;
                }
                idim += 1;
            }

            // 3D subspaces (triples of coordinates).
            for ii in 1..=num3dss {
                if !is_ci {
                    break;
                }
                let idx = idim;
                let eps = local_tolerance[idx - 1] * 0.01;
                let o1 = base1 + (num1dss + 2 * num2dss) as usize + (ii as usize - 1) * 3;
                let o2 = half + base2 + (num1dss + 2 * num2dss) as usize + (ii as usize - 1) * 3;
                let v1x = result[o1] * facteur1;
                let v1y = result[o1 + 1] * facteur1;
                let v1z = result[o1 + 2] * facteur1;
                let v2x = result[o2] * facteur2;
                let v2y = result[o2 + 1] * facteur2;
                let v2z = result[o2 + 2] * facteur2;
                let diff = (v1x - v2x).abs() + (v1y - v2y).abs() + (v1z - v2z).abs();
                let moy = (v1x + v2x).abs() + (v1y + v2y).abs() + (v1z + v2z).abs();
                if diff > moy * 1.0e-9 {
                    prec[idim] = diff * normal1;
                    suivant[idim] = diff * normal2;
                    if prec[idim] > eps || suivant[idim] > eps {
                        is_ci = false;
                    }
                } else {
                    prec[idim] = 0.0;
                    suivant[idim] = 0.0;
                }
                idim += 1;
            }

            // If it's good, update everything.
            if is_ci {
                continuity[icurve as usize] = iordre;
                let index = ((icurve - 1) * nbspace as i32) as usize;
                for ii in 0..nbspace {
                    error_max[index + ii] += prec[ii + 1];
                    error_max[index + nbspace + ii] += suivant[ii + 1];
                }
            }
        }
    }
}

/// OCCT AdvApprox_ApproxAFunction::Approximation (static, L364-599) —
/// approximates a non-polynomial function by polynomial curves with interval
/// cutting.  Outputs are the flat arrays described by the OCCT signature.
#[allow(clippy::too_many_arguments)]
pub fn approximation(
    total_dimension: usize,
    total_num_ss: usize,
    local_dimension: &[i32],
    first: f64,
    last: f64,
    evaluator: &mut dyn EvaluatorFunction,
    cut_tool: &dyn Cutting,
    continuity_order: i32,
    num_max_coeffs: usize,
    max_segments: usize,
    local_tolerances_array: &[f64],
    code_precis: i32,
    num_curves: &mut i32,
    num_coeff_per_curve_array: &mut [i32],
    coefficient_array: &mut [f64],
    intervals_array: &mut [f64],
    error_max_array: &mut [f64],
    average_error_array: &mut [f64],
    error_code: &mut i32,
) {
    let mut is_cut = false;

    *error_code = 0;
    for v in coefficient_array.iter_mut() {
        *v = 0.0;
    }

    //-------------------------- Input validation ------------------
    if max_segments < 1 || (last - first).abs() < 1.0e-9 {
        *error_code = 1;
        return;
    }

    //--> The total dimension must be the sum of the subspace dimensions.
    let mut idim = 0usize;
    for i in 0..total_num_ss {
        idim += local_dimension[i] as usize;
    }
    if idim != total_dimension {
        *error_code = 1;
        return;
    }
    let continuity = match continuity_order {
        0 => GeomAbsShape::C0,
        1 => GeomAbsShape::C1,
        2 => GeomAbsShape::C2,
        _ => panic!("Standard_ConstructionError"),
    };

    //--------------------- Choice of number of points ----------------------
    let mut nb_gauss_points = 0usize;
    let mut work_degree = 0usize;
    jacobi_parameters(
        continuity,
        num_max_coeffs - 1,
        code_precis,
        &mut nb_gauss_points,
        &mut work_degree,
    );

    //------------------ Initialization of cutting management ---------
    // TABINT = IntervalsArray: [0] = First, [1] = Last, growing.
    intervals_array[0] = first;
    intervals_array[1] = last;
    let mut nupil = 1usize;
    *num_curves = 0;

    // ********************************************************************
    //                      APPROXIMATION WITH CUTTING
    // ********************************************************************
    let jacobi_base = JacobiPolynomial::new(work_degree, continuity);
    let max_degree = num_max_coeffs - 1;
    {
        let mut approx = SimpleApprox::new(
            total_dimension,
            total_num_ss,
            continuity,
            work_degree,
            nb_gauss_points,
            &jacobi_base,
        );

        // TabINT stack; slots [num_curves ..= nupil].
        while nupil - (*num_curves as usize) != 0 {
            //---- Compute the approximation curve in the Jacobi basis -----
            approx.perform(
                evaluator,
                local_dimension,
                local_tolerances_array,
                intervals_array[*num_curves as usize],
                intervals_array[*num_curves as usize + 1],
                max_degree,
            );
            if !approx.is_done() {
                *error_code = 1;
                return;
            }

            //---------- Compute the curve degree and max error ----------
            // (NumCoeffPerCurveArray(NumCurves+1) = 0 in OCCT.)

            //    The error must be satisfied on all subspaces, otherwise split.
            let mut max_max_err = true;
            for is in 0..total_num_ss {
                if approx.max_error(is + 1) > local_tolerances_array[is] {
                    max_max_err = false;
                    break;
                }
            }

            if max_max_err {
                *num_curves += 1;
            } else {
                //-> ...otherwise try to split the current interval in 2...
                let large = cut_tool.value(
                    intervals_array[*num_curves as usize],
                    intervals_array[*num_curves as usize + 1],
                );

                if nupil < max_segments && large.is_some() {
                    let tmil = large.unwrap();
                    is_cut = true; // Now we know it!
                    // Shift right: IDEB1[iI] = IDEB[iI] for iI = ILONG..0,
                    // then IDEB[0] = TMIL.
                    let ilong = nupil - *num_curves as usize - 1;
                    let ideb = *num_curves as usize + 1;
                    for ii in (0..=ilong).rev() {
                        intervals_array[ideb + 1 + ii] = intervals_array[ideb + ii];
                    }
                    intervals_array[ideb] = tmil;
                    nupil += 1;
                    //--> ... and start over.
                    continue;
                } else {
                    //--> If the stack is full...
                    *num_curves += 1;
                }
            }

            for is in 0..total_num_ss {
                error_max_array[is + (*num_curves as usize - 1) * total_num_ss] =
                    approx.max_error(is + 1);
                average_error_array[is + (*num_curves as usize - 1) * total_num_ss] =
                    approx.average_error(is + 1);
            }

            let jac_coeff: Vec<f64> = approx.coefficients().to_vec();
            let mut the_deg = approx.degree() as i32;
            if is_cut && (the_deg < 2 * continuity_order + 1) {
                // To avoid noisy derivatives at the ends, and maintain
                // correct continuity on the resulting BSpline.
                the_deg = 2 * continuity_order + 1;
            }

            num_coeff_per_curve_array[(*num_curves - 1) as usize] = the_deg + 1;
            let mut coefficients = vec![0.0f64; ((the_deg + 1) * total_dimension as i32) as usize];
            jacobi_base.to_coefficients(
                total_dimension,
                the_deg as usize,
                &jac_coeff,
                &mut coefficients,
            );
            let f = (the_deg + 1) as usize * total_dimension;
            let mut j = ((*num_curves as i32 - 1) * total_dimension as i32 * num_max_coeffs as i32) as usize;
            for i in 0..f {
                coefficient_array[j] = coefficients[i];
                j += 1;
            }
        }
    }
}

/// OCCT AdvApprox_ApproxAFunction.
pub struct ApproxAFunction {
    num_sub_spaces: [i32; 3],
    first: f64,
    last: f64,
    continuity: GeomAbsShape,
    max_degree: i32,
    max_segments: i32,
    done: bool,
    has_result: bool,
    knots: Vec<f64>,
    mults: Vec<i32>,
    degree: i32,
    poles3d: Vec<Vec<f64>>, // [pole_idx][ss_idx] flattened coordinate triples? kept as OCCT Array2 (row-major poles x dimension)
    poles3d_ncols: usize,
    max_error_3d: Vec<f64>,
    average_error_3d: Vec<f64>,
}

impl ApproxAFunction {
    /// OCCT AdvApprox_ApproxAFunction ctor with the default dichotomy
    /// cutting, followed by Perform.  Only the 3D-subspace configuration used
    /// by the helix pipeline is stored (Num1DSS = 0, Num2DSS = 0).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        num1dss: i32,
        num2dss: i32,
        num3dss: i32,
        one_d_tol: Option<&[f64]>,
        two_d_tol: Option<&[f64]>,
        three_d_tol: Option<&[f64]>,
        first: f64,
        last: f64,
        continuity: GeomAbsShape,
        max_deg: i32,
        max_seg: i32,
        func: &mut dyn EvaluatorFunction,
    ) -> Self {
        let mut this = ApproxAFunction {
            num_sub_spaces: [num1dss, num2dss, num3dss],
            first,
            last,
            continuity,
            max_degree: max_deg,
            max_segments: max_seg,
            done: false,
            has_result: false,
            knots: Vec::new(),
            mults: Vec::new(),
            degree: 0,
            poles3d: Vec::new(),
            poles3d_ncols: 0,
            max_error_3d: Vec::new(),
            average_error_3d: Vec::new(),
        };
        let mut cut = super::cutting::DichoCutting;
        this.perform(
            num1dss,
            num2dss,
            num3dss,
            one_d_tol,
            two_d_tol,
            three_d_tol,
            func,
            &cut,
        );
        this
    }

    /// OCCT AdvApprox_ApproxAFunction::Perform.
    #[allow(clippy::too_many_arguments)]
    fn perform(
        &mut self,
        num1dss: i32,
        num2dss: i32,
        num3dss: i32,
        one_d_tol: Option<&[f64]>,
        two_d_tol: Option<&[f64]>,
        three_d_tol: Option<&[f64]>,
        func: &mut dyn EvaluatorFunction,
        cut_tool: &dyn Cutting,
    ) {
        if num1dss < 0 || num2dss < 0 || num3dss < 0 || num1dss + num2dss + num3dss <= 0
            || self.last < self.first || self.max_degree < 1 || self.max_segments < 0
        {
            panic!("Standard_ConstructionError");
        }
        if self.max_degree > 14 {
            self.max_degree = 14;
        }

        self.num_sub_spaces = [num1dss, num2dss, num3dss];
        let total_num_ss = (num1dss + num2dss + num3dss) as usize;
        let total_dimension = (self.num_sub_spaces[0] + 2 * self.num_sub_spaces[1] + 3 * self.num_sub_spaces[2]) as usize;

        let mut continuity_order = 0i32;
        match self.continuity {
            GeomAbsShape::C0 => continuity_order = 0,
            GeomAbsShape::C1 => continuity_order = 1,
            GeomAbsShape::C2 => continuity_order = 2,
            _ => panic!("Standard_ConstructionError"),
        }
        let approx_start_end = [-1.0f64, 1.0f64];
        let num_max_coeffs = ((self.max_degree + 1) as usize).max((2 * continuity_order + 2) as usize);
        self.max_degree = num_max_coeffs as i32 - 1;
        let code_precis = 1;

        // LocalDimension / LocalTolerances.
        let mut local_dimension = vec![0i32; total_num_ss];
        let mut local_tolerances = vec![0.0f64; total_num_ss];
        let mut index = 0usize;
        for jj in 0..self.num_sub_spaces[0] as usize {
            local_tolerances[index] = one_d_tol.unwrap()[jj];
            local_dimension[index] = 1;
            index += 1;
        }
        for jj in 0..self.num_sub_spaces[1] as usize {
            local_tolerances[index] = two_d_tol.unwrap()[jj];
            local_dimension[index] = 2;
            index += 1;
        }
        for jj in 0..self.num_sub_spaces[2] as usize {
            local_tolerances[index] = three_d_tol.unwrap()[jj];
            local_dimension[index] = 3;
            index += 1;
        }

        // Output.
        let mut error_code = 0i32;
        let mut num_curves = 0i32;
        let size = self.max_segments as usize * num_max_coeffs * total_dimension;
        let mut num_coeff_per_curve = vec![0i32; self.max_segments as usize];
        let mut local_coefficients = vec![0.0f64; size];
        let mut intervals = vec![0.0f64; self.max_segments as usize + 1];
        let mut error_max = vec![0.0f64; self.max_segments as usize * total_num_ss];
        let mut average_error = vec![0.0f64; self.max_segments as usize * total_num_ss];

        approximation(
            total_dimension,
            total_num_ss,
            &local_dimension,
            self.first,
            self.last,
            func,
            cut_tool,
            continuity_order,
            num_max_coeffs,
            self.max_segments as usize,
            &local_tolerances,
            code_precis,
            &mut num_curves,
            &mut num_coeff_per_curve,
            &mut local_coefficients,
            &mut intervals,
            &mut error_max,
            &mut average_error,
            &mut error_code,
        );

        if error_code == 0 || error_code == -1 {
            // Everything OK, or a result with one error above the tolerance.
            let mut tab_continuity = vec![0i32; num_curves as usize];
            let mut polynomial_intervals = vec![0.0f64; num_curves as usize * 2];
            for ii in 0..num_curves as usize {
                // Force a minimum degree of 1 (PRO5474).
                num_coeff_per_curve[ii] = num_coeff_per_curve[ii].max(2);
                polynomial_intervals[ii * 2] = approx_start_end[0];
                polynomial_intervals[ii * 2 + 1] = approx_start_end[1];
            }

            prepare_convert(
                num_curves,
                self.max_degree,
                continuity_order,
                num1dss,
                num2dss,
                num3dss,
                &num_coeff_per_curve,
                &mut local_coefficients,
                &polynomial_intervals,
                &intervals,
                &local_tolerances,
                &mut error_max,
                &mut tab_continuity,
            );

            let a_converter = ConvertCompPolynomialToPoles::from_arrays(
                num_curves as usize,
                total_dimension,
                self.max_degree as usize,
                &tab_continuity,
                &num_coeff_per_curve,
                &local_coefficients,
                &polynomial_intervals,
                &intervals,
            );

            if a_converter.is_done() {
                let poles = a_converter.poles(); // row-major (nb_poles x dimension)
                let nb_poles_rows = a_converter.nb_poles();
                self.knots = a_converter.knots_vec();
                self.mults = a_converter.multiplicities_vec();
                self.degree = a_converter.degree() as i32;
                let dim_index = 0usize; // no 1D/2D subspaces stored
                if self.num_sub_spaces[2] > 0 {
                    let nss = self.num_sub_spaces[2] as usize;
                    self.poles3d = vec![vec![0.0f64; nss * 3]; nb_poles_rows];
                    self.poles3d_ncols = nss * 3;
                    self.max_error_3d = vec![0.0f64; nss];
                    self.average_error_3d = vec![0.0f64; nss];
                    for ii in 0..nb_poles_rows {
                        for jj in 0..nss {
                            let local_index = dim_index + jj * 3;
                            for kk in 0..3 {
                                self.poles3d[ii][jj * 3 + kk] =
                                    poles[ii * total_dimension + local_index + kk];
                            }
                        }
                    }
                    for jj in 0..nss {
                        let mut error_value = 0.0f64;
                        for ii in 0..num_curves as usize {
                            let local_index = ii * total_num_ss;
                            error_value = error_max[local_index + jj].max(error_value);
                        }
                        self.max_error_3d[jj] = error_value;
                    }
                    for jj in 0..nss {
                        let mut error_value = 0.0f64;
                        for ii in 0..num_curves as usize {
                            let local_index = ii * total_num_ss;
                            error_value += average_error[local_index + jj];
                        }
                        error_value /= num_curves as f64;
                        self.average_error_3d[jj] = error_value;
                    }
                }
                if error_code == 0 {
                    self.done = true;
                    self.has_result = true;
                } else if error_code == -1 {
                    self.has_result = true;
                }
            }
        }
    }

    /// OCCT AdvApprox_ApproxAFunction::IsDone.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT AdvApprox_ApproxAFunction::HasResult.
    pub fn has_result(&self) -> bool {
        self.has_result
    }

    /// OCCT AdvApprox_ApproxAFunction::NbPoles.
    pub fn nb_poles(&self) -> usize {
        if self.done || self.has_result {
            return nb_poles(self.degree as usize, false, &self.mults);
        }
        0
    }

    /// The 3D poles of subspace `index` (1-based) as a flat
    /// [pole][coord] array (row-major), matching OCCT Poles(Index, P).
    pub fn poles_flat(&self, index: usize) -> Vec<f64> {
        let mut out = vec![0.0f64; self.poles3d.len() * 3];
        for (ii, row) in self.poles3d.iter().enumerate() {
            let o = (index - 1) * 3;
            out[ii * 3..ii * 3 + 3].copy_from_slice(&row[o..o + 3]);
        }
        out
    }

    /// OCCT AdvApprox_ApproxAFunction::Degree.
    pub fn degree(&self) -> i32 {
        self.degree
    }

    /// OCCT AdvApprox_ApproxAFunction::Knots.
    pub fn knots_vec(&self) -> &[f64] {
        &self.knots
    }

    /// OCCT AdvApprox_ApproxAFunction::Multiplicities.
    pub fn multiplicities_vec(&self) -> &[i32] {
        &self.mults
    }

    /// OCCT AdvApprox_ApproxAFunction::MaxError(D, Index) — only the 3D
    /// subspace is populated (D == 3).
    pub fn max_error_at(&self, d: usize, index: usize) -> f64 {
        assert!(d == 3, "AdvApprox: only 3D subspace errors are stored");
        self.max_error_3d[index - 1]
    }

    /// OCCT AdvApprox_ApproxAFunction::AverageError(D, Index) — 3D subspace.
    pub fn average_error_at(&self, d: usize, index: usize) -> f64 {
        assert!(d == 3, "AdvApprox: only 3D subspace errors are stored");
        self.average_error_3d[index - 1]
    }

    /// The poles matrix row-major as (nb_poles x dimension) rows of xyz
    /// triples per subspace — exposed for HelixGeom_Tools.
    pub fn poles3d_rows(&self) -> &Vec<Vec<f64>> {
        &self.poles3d
    }

    /// Poles column count (nss * 3).
    pub fn poles3d_ncols(&self) -> usize {
        self.poles3d_ncols
    }

    /// Number of subspace rows in the poles matrix.
    pub fn poles3d_nrows(&self) -> usize {
        self.poles3d.len()
    }
}
