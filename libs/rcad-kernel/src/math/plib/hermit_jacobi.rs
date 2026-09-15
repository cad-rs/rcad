//! OCCT PLib_HermitJacobi (FoundationClasses/TKMath/PLib).
//!
//! 1:1 translation of `PLib_HermitJacobi.hxx` (L17-142) and
//! `PLib_HermitJacobi.cxx` (L1-373, including the D0123 basis assembly,
//! FORTRAN MPOBAS port).  There is no `PLib_Base` abstract class in this
//! OCCT version - `PLib_HermitJacobi` is a plain class composing a
//! `PLib_JacobiPolynomial` member (hxx L139), so the translation mirrors
//! that composition directly.
//!
//! Container mapping: OCCT `NCollection_Array1<double>` basis arrays are
//! `&mut [f64]` slices (all OCCT call sites pass 0-based arrays);
//! `math_Matrix` locals use the bound-aware `math_matrix::Matrix` mirror;
//! `NCollection_LocalArray` stack buffers become fixed-size arrays.

use std::sync::OnceLock;

use super::super::math_matrix::Matrix;
use super::super::plib::{hermite_coefficients, no_derivative_eval_polynomial_flat};
use super::super::GeomAbsShape;
use super::jacobi_polynomial::JacobiPolynomial;

/// W coefficients for C0 continuity (NivConstr = 0, DegreeH = 1)
/// W(t) = (1 - t^2) (PLib_HermitJacobi.cxx L28).
const W_COEFF_C0: [f64; 3] = [1.0, 0.0, -1.0];

/// W coefficients for C1 continuity (NivConstr = 1, DegreeH = 3)
/// W(t) = (1 - t^2)^2 = 1 - 2t^2 + t^4 (PLib_HermitJacobi.cxx L32).
const W_COEFF_C1: [f64; 5] = [1.0, 0.0, -2.0, 0.0, 1.0];

/// W coefficients for C2 continuity (NivConstr = 2, DegreeH = 5)
/// W(t) = (1 - t^2)^3 = 1 - 3t^2 + 3t^4 - t^6 (PLib_HermitJacobi.cxx L36).
const W_COEFF_C2: [f64; 7] = [1.0, 0.0, -3.0, 0.0, 3.0, 0.0, -1.0];

/// OCCT GetWCoefficients (PLib_HermitJacobi.cxx L38-51): returns the W
/// polynomial coefficient table for the given NivConstr.
#[inline]
fn get_w_coefficients(the_niv_constr: i32) -> &'static [f64] {
    match the_niv_constr {
        0 => &W_COEFF_C0,
        1 => &W_COEFF_C1,
        2 => &W_COEFF_C2,
        // Fallback, should never happen
        _ => &W_COEFF_C0,
    }
}

/// OCCT GetHermiteMatrix_C0 (PLib_HermitJacobi.cxx L53-61): the C++ static
/// local becomes a OnceLock-initialized matrix.
fn get_hermite_matrix_c0() -> &'static Matrix {
    static A_MATRIX: OnceLock<Matrix> = OnceLock::new();
    A_MATRIX.get_or_init(|| {
        let mut a_result = Matrix::new(1, 2, 1, 2);
        // OCCT discards the bool result of PLib::HermiteCoefficients here.
        let _ = hermite_coefficients(-1.0, 1.0, 0, 0, &mut a_result);
        a_result
    })
}

/// OCCT GetHermiteMatrix_C1 (PLib_HermitJacobi.cxx L63-70).
fn get_hermite_matrix_c1() -> &'static Matrix {
    static A_MATRIX: OnceLock<Matrix> = OnceLock::new();
    A_MATRIX.get_or_init(|| {
        let mut a_result = Matrix::new(1, 4, 1, 4);
        let _ = hermite_coefficients(-1.0, 1.0, 1, 1, &mut a_result);
        a_result
    })
}

/// OCCT GetHermiteMatrix_C2 (PLib_HermitJacobi.cxx L72-79).
fn get_hermite_matrix_c2() -> &'static Matrix {
    static A_MATRIX: OnceLock<Matrix> = OnceLock::new();
    A_MATRIX.get_or_init(|| {
        let mut a_result = Matrix::new(1, 6, 1, 6);
        let _ = hermite_coefficients(-1.0, 1.0, 2, 2, &mut a_result);
        a_result
    })
}

/// OCCT GetHermiteMatrix (PLib_HermitJacobi.cxx L83-96).
#[inline]
fn get_hermite_matrix(the_niv_constr: i32) -> &'static Matrix {
    match the_niv_constr {
        0 => get_hermite_matrix_c0(),
        1 => get_hermite_matrix_c1(),
        2 => get_hermite_matrix_c2(),
        // Fallback, should never happen
        _ => get_hermite_matrix_c0(),
    }
}

/// OCCT PLib::EvalPolynomial (PLib.cxx L945-1016, with eval_poly1 L746-773,
/// eval_poly2 L779-812 and the general default branch L1000-1014) - value
/// and derivatives of a polynomial with ASCENDING canonical coefficients
/// ([0] = X^0 coefficient ... [d * Dimension] = X^Degree coefficient).
///
/// Hosting note: the kernel's existing `math::plib::eval_polynomial_flat`
/// walks the slice from the LOW index as if it held the highest-degree
/// coefficient (descending storage), so it cannot be reused for the
/// ascending Hermite-matrix rows consumed here; this is the faithful
/// ascending-storage form.
fn eval_polynomial(
    par: f64,
    derivative_request: i32,
    degree: i32,
    dimension: i32,
    polynomial_coeff: &[f64],
    results: &mut [f64],
) {
    let dim = dimension as usize;
    // OCCT: aCoeffs points at the highest-degree coefficient block.
    let mut a_coeffs = degree as usize * dim;

    match derivative_request {
        1 => {
            // eval_poly1
            let mut local0 = vec![0.0f64; dim];
            let mut local1 = vec![0.0f64; dim];

            for i in 0..dim {
                local0[i] = polynomial_coeff[a_coeffs + i];
                local1[i] = 0.0;
            }

            for _a_deg in 0..degree {
                a_coeffs -= dim;
                for i in 0..dim {
                    let a_val = local0[i];
                    local1[i] = local1[i] * par + a_val;
                    local0[i] = a_val * par + polynomial_coeff[a_coeffs + i];
                }
            }

            for i in 0..dim {
                results[i] = local0[i];
                results[dim + i] = local1[i];
            }
        }
        2 => {
            // eval_poly2
            let mut local0 = vec![0.0f64; dim];
            let mut local1 = vec![0.0f64; dim];
            let mut local2 = vec![0.0f64; dim];

            for i in 0..dim {
                local0[i] = polynomial_coeff[a_coeffs + i];
                local1[i] = 0.0;
                local2[i] = 0.0;
            }

            for _a_deg in 0..degree {
                a_coeffs -= dim;
                for i in 0..dim {
                    let a_d1 = local1[i];
                    let a_val = local0[i];
                    local2[i] = local2[i] * par + a_d1 * 2.0;
                    local1[i] = a_d1 * par + a_val;
                    local0[i] = a_val * par + polynomial_coeff[a_coeffs + i];
                }
            }

            for i in 0..dim {
                results[i] = local0[i];
                results[dim + i] = local1[i];
                results[2 * dim + i] = local2[i];
            }
        }
        _ => {
            // General case (DerivativeRequest == 0 and DerivativeRequest > 2)
            let res_size = (1 + derivative_request) as usize * dim;
            for v in results[..res_size].iter_mut() {
                *v = 0.0;
            }

            for _a_deg in 0..=degree {
                let mut a_ptr = res_size - dim;
                // Calculating derivatives of the polynomial
                for a_deriv in (1..=derivative_request).rev() {
                    let an_original = a_ptr - dim;
                    for ind in 0..dim {
                        results[a_ptr + ind] = results[a_ptr + ind] * par
                            + results[an_original + ind] * a_deriv as f64;
                    }
                    a_ptr = an_original;
                }
                // Calculating the value of the polynomial
                for ind in 0..dim {
                    results[a_ptr + ind] =
                        results[a_ptr + ind] * par + polynomial_coeff[a_coeffs + ind];
                }
                if a_coeffs >= dim {
                    a_coeffs -= dim;
                }
            }
        }
    }
}

/// OCCT PLib_HermitJacobi - Hermit-Jacobi basis relative to an order of
/// constraint (PLib_HermitJacobi.hxx L27-57):
/// P(t) = H(t) + W(t) * Q(t), W(t) = (1-t^2)^(2*nivConstr+2).
#[derive(Debug, Clone)]
pub struct HermitJacobi {
    my_jacobi: JacobiPolynomial,
}

impl HermitJacobi {
    /// OCCT PLib_HermitJacobi::PLib_HermitJacobi (PLib_HermitJacobi.cxx
    /// L101-104).
    pub fn new(the_work_degree: i32, the_constraint_order: GeomAbsShape) -> Self {
        HermitJacobi {
            my_jacobi: JacobiPolynomial::new(the_work_degree, the_constraint_order),
        }
    }

    /// OCCT PLib_HermitJacobi::WorkDegree (hxx L123).
    #[inline]
    pub fn work_degree(&self) -> i32 {
        self.my_jacobi.work_degree()
    }

    /// OCCT PLib_HermitJacobi::NivConstr (hxx L126).
    #[inline]
    pub fn niv_constr(&self) -> i32 {
        self.my_jacobi.niv_constr()
    }

    /// OCCT PLib_HermitJacobi::MaxError (PLib_HermitJacobi.cxx L108-113).
    pub fn max_error(&self, dimension: i32, herm_jac_coeff: &[f64], new_degree: i32) -> f64 {
        self.my_jacobi.max_error(dimension, herm_jac_coeff, new_degree)
    }

    /// OCCT PLib_HermitJacobi::ReduceDegree (PLib_HermitJacobi.cxx L117-125).
    pub fn reduce_degree(
        &self,
        dimension: i32,
        max_degree: i32,
        tol: f64,
        herm_jac_coeff: &[f64],
        new_degree: &mut i32,
        max_error: &mut f64,
    ) {
        self.my_jacobi
            .reduce_degree(dimension, max_degree, tol, herm_jac_coeff, new_degree, max_error);
    }

    /// OCCT PLib_HermitJacobi::AverageError (PLib_HermitJacobi.cxx L129-134).
    pub fn average_error(&self, dimension: i32, herm_jac_coeff: &[f64], new_degree: i32) -> f64 {
        self.my_jacobi.average_error(dimension, herm_jac_coeff, new_degree)
    }

    /// OCCT PLib_HermitJacobi::ToCoefficients (PLib_HermitJacobi.cxx
    /// L138-188).  Converts the polynomial P(t) = H(t) + W(t) Q(t) in the
    /// canonical base.  Both coefficient arrays use the Lower() == 0 slice
    /// convention.
    pub fn to_coefficients(
        &self,
        dimension: i32,
        degree: i32,
        herm_jac_coeff: &[f64],
        coefficients: &mut [f64],
    ) {
        let a_niv_constr = self.niv_constr();
        let a_degree_h = 2 * a_niv_constr + 1;
        let ibeg_hjc = 0usize; // OCCT: HermJacCoeff.Lower() (== 0 convention)

        let a_hermite_matrix = get_hermite_matrix(a_niv_constr);

        let mut aux_coeff = vec![0.0f64; ((degree + 1) * dimension) as usize];

        for k in 0..=a_degree_h {
            let kdim = (k * dimension) as usize;
            for i in 0..=a_niv_constr {
                let h1 = a_hermite_matrix.get(i + 1, k + 1);
                let h2 = a_hermite_matrix.get(i + a_niv_constr + 2, k + 1);
                let i1 = ibeg_hjc + (i * dimension) as usize;
                let i2 = ibeg_hjc + ((i + a_niv_constr + 1) * dimension) as usize;

                for idim in 0..dimension {
                    let idimu = idim as usize;
                    aux_coeff[idimu + kdim] +=
                        herm_jac_coeff[i1 + idimu] * h1 + herm_jac_coeff[i2 + idimu] * h2;
                }
            }
        }
        let kdim = ((degree + 1) * dimension) as usize;
        for k in ((a_degree_h + 1) * dimension) as usize..kdim {
            aux_coeff[k] = herm_jac_coeff[ibeg_hjc + k];
        }

        if degree > a_degree_h {
            self.my_jacobi.to_coefficients(dimension, degree, &aux_coeff, coefficients);
        } else {
            for (c, a) in coefficients.iter_mut().zip(aux_coeff.iter()) {
                *c = *a;
            }
        }
    }

    /// OCCT PLib_HermitJacobi::D0123 (PLib_HermitJacobi.cxx L195-335),
    /// common part of D0, D1, D2, D3 (FORTRAN MPOBAS port).
    ///
    /// Architecture deviation: OCCT passes the SAME array several times
    /// (`D0: D0123(0, U, V, V, V, V)`), which the C++ aliasing rules allow
    /// and Rust `&mut` forbids.  The derivative outputs therefore arrive as
    /// `Option`, consumed only for the derivative order actually requested;
    /// the public D0..D3 wrappers below reproduce the OCCT call forms.
    fn d0123(
        &self,
        n_deriv: i32,
        u: f64,
        basis_value: &mut [f64],
        basis_d1: Option<&mut [f64]>,
        basis_d2: Option<&mut [f64]>,
        basis_d3: Option<&mut [f64]>,
    ) {
        // NCollection_LocalArray<double> jac0(4 * 20) etc.
        let mut jac0 = [0.0f64; 4 * 20];
        let mut jac1 = [0.0f64; 4 * 20];
        let mut jac2 = [0.0f64; 4 * 20];
        let mut jac3 = [0.0f64; 4 * 20];
        // NCollection_LocalArray<double> wvalues(4); WValues(0, NDeriv).
        let mut w_values = [0.0f64; 4];

        let a_niv_constr = self.niv_constr();
        let a_work_degree = self.work_degree();
        let a_degree_h = 2 * a_niv_constr + 1;
        // ibeg0 = ibeg1 = ibeg2 = ibeg3 = 0 (Lower() == 0 slice convention).
        let a_jac_degree = a_work_degree - a_degree_h - 1;

        // Consume the derivative outputs for the requested order (see the
        // Option-based signature note above); the dummy fallbacks are never
        // dereferenced because every use below is guarded by NDeriv.
        let mut empty1: Vec<f64> = Vec::new();
        let mut empty2: Vec<f64> = Vec::new();
        let mut empty3: Vec<f64> = Vec::new();
        let mut basis_d1 = basis_d1;
        let mut basis_d2 = basis_d2;
        let mut basis_d3 = basis_d3;
        let basis_d1: &mut [f64] = if n_deriv >= 1 {
            basis_d1.take().expect("BasisD1 required for NDeriv >= 1")
        } else {
            &mut empty1[..]
        };
        let basis_d2: &mut [f64] = if n_deriv >= 2 {
            basis_d2.take().expect("BasisD2 required for NDeriv >= 2")
        } else {
            &mut empty2[..]
        };
        let basis_d3: &mut [f64] = if n_deriv >= 3 {
            basis_d3.take().expect("BasisD3 required for NDeriv == 3")
        } else {
            &mut empty3[..]
        };

        let a_hermite_matrix = get_hermite_matrix(a_niv_constr);

        let a_jac_len = (a_jac_degree.max(0) + 1) as usize;
        let jac_value0: &mut [f64] = &mut jac0[..a_jac_len];
        let w_values_slice: &mut [f64] = &mut w_values[..(n_deriv + 1) as usize];

        // Evaluation des polynomes d'hermite.
        // HermitValues(0, aDegreeH, 0, NDeriv); logical row i (lower 0) is
        // storage row i; the Hermite matrix logical row (i + 1) is storage
        // row i.
        let mut hermit_values = Matrix::new(0, a_degree_h, 0, n_deriv);
        if n_deriv == 0 {
            for i in 0..=a_degree_h {
                no_derivative_eval_polynomial_flat(
                    u,
                    a_degree_h,
                    1,
                    a_degree_h,
                    &a_hermite_matrix.data.m[i as usize],
                    &mut hermit_values.data.m[i as usize],
                );
            }
        } else {
            for i in 0..=a_degree_h {
                eval_polynomial(
                    u,
                    n_deriv,
                    a_degree_h,
                    1,
                    &a_hermite_matrix.data.m[i as usize],
                    &mut hermit_values.data.m[i as usize],
                );
            }
        }

        // Evaluation des polynomes de Jaccobi
        if a_jac_degree >= 0 {
            match n_deriv {
                0 => {
                    self.my_jacobi.d0(u, jac_value0);
                }
                1 => {
                    self.my_jacobi.d1(u, jac_value0, &mut jac1[..a_jac_len]);
                }
                2 => {
                    self.my_jacobi
                        .d2(u, jac_value0, &mut jac1[..a_jac_len], &mut jac2[..a_jac_len]);
                }
                _ => {
                    self.my_jacobi.d3(
                        u,
                        jac_value0,
                        &mut jac1[..a_jac_len],
                        &mut jac2[..a_jac_len],
                        &mut jac3[..a_jac_len],
                    );
                }
            }

            // Evaluation de W(t)
            let a_w_coeff = get_w_coefficients(a_niv_constr);
            if n_deriv == 0 {
                no_derivative_eval_polynomial_flat(
                    u,
                    a_degree_h + 1,
                    1,
                    a_degree_h + 1,
                    a_w_coeff,
                    w_values_slice,
                );
            } else {
                eval_polynomial(u, n_deriv, a_degree_h + 1, 1, a_w_coeff, w_values_slice);
            }
        }

        // Evaluation a l'ordre 0
        for i in 0..=a_degree_h {
            basis_value[i as usize] = hermit_values.get(i, 0);
        }
        let w0 = w_values[0];
        for i in (a_degree_h + 1)..=a_work_degree {
            let j = (i - a_degree_h - 1) as usize;
            basis_value[i as usize] = w0 * jac0[j];
        }

        // Evaluation a l'ordre 1
        if n_deriv >= 1 {
            let w1 = w_values[1];
            for i in 0..=a_degree_h {
                basis_d1[i as usize] = hermit_values.get(i, 1);
            }
            for i in (a_degree_h + 1)..=a_work_degree {
                let j = (i - a_degree_h - 1) as usize;
                basis_d1[i as usize] = w0 * jac1[j] + w1 * jac0[j];
            }
            // Evaluation a l'ordre 2
            if n_deriv >= 2 {
                let w2 = w_values[2];
                for i in 0..=a_degree_h {
                    basis_d2[i as usize] = hermit_values.get(i, 2);
                }
                for i in (a_degree_h + 1)..=a_work_degree {
                    let j = (i - a_degree_h - 1) as usize;
                    basis_d2[i as usize] = w0 * jac2[j] + 2.0 * w1 * jac1[j] + w2 * jac0[j];
                }

                // Evaluation a l'ordre 3
                if n_deriv == 3 {
                    let w3 = w_values[3];
                    for i in 0..=a_degree_h {
                        basis_d3[i as usize] = hermit_values.get(i, 3);
                    }
                    for i in (a_degree_h + 1)..=a_work_degree {
                        let j = (i - a_degree_h - 1) as usize;
                        basis_d3[i as usize] =
                            w0 * jac3[j] + w3 * jac0[j] + 3.0 * (w1 * jac2[j] + w2 * jac1[j]);
                    }
                }
            }
        }
    }

    /// OCCT PLib_HermitJacobi::D0 (PLib_HermitJacobi.cxx L339-342).
    pub fn d0(&self, u: f64, basis_value: &mut [f64]) {
        self.d0123(0, u, basis_value, None, None, None);
    }

    /// OCCT PLib_HermitJacobi::D1 (PLib_HermitJacobi.cxx L346-351).
    pub fn d1(&self, u: f64, basis_value: &mut [f64], basis_d1: &mut [f64]) {
        self.d0123(1, u, basis_value, Some(basis_d1), None, None);
    }

    /// OCCT PLib_HermitJacobi::D2 (PLib_HermitJacobi.cxx L355-361).
    pub fn d2(&self, u: f64, basis_value: &mut [f64], basis_d1: &mut [f64], basis_d2: &mut [f64]) {
        self.d0123(2, u, basis_value, Some(basis_d1), Some(basis_d2), None);
    }

    /// OCCT PLib_HermitJacobi::D3 (PLib_HermitJacobi.cxx L365-372).
    pub fn d3(
        &self,
        u: f64,
        basis_value: &mut [f64],
        basis_d1: &mut [f64],
        basis_d2: &mut [f64],
        basis_d3: &mut [f64],
    ) {
        self.d0123(
            3,
            u,
            basis_value,
            Some(basis_d1),
            Some(basis_d2),
            Some(basis_d3),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT PLib_HermitJacobi D0 basis, hand-derived for NivConstr = 0
    /// (DegreeH = 1, W(t) = 1 - t^2), WorkDegree = 2 (Jacobi degree 0):
    ///   Basis[0] = H00(u) = (1 - u)/2      (value 1 at u = -1)
    ///   Basis[1] = H01(u) = (1 + u)/2      (value 1 at u = +1)
    ///   Basis[2] = W(u) * J0(u) = (1 - u^2) * sqrt(15/16)
    /// at u = 1/2: [0.25, 0.75, 0.75*sqrt(15/16)].
    #[test]
    fn hermit_jacobi_d0_hand() {
        let base = HermitJacobi::new(2, GeomAbsShape::C0);

        let mut v = [0.0f64; 3];
        base.d0(0.5, &mut v);

        assert!((v[0] - 0.25).abs() < 1e-14, "v0 = {}", v[0]);
        assert!((v[1] - 0.75).abs() < 1e-14, "v1 = {}", v[1]);
        let expected2 = 0.75 * (15.0f64 / 16.0).sqrt();
        assert!((v[2] - expected2).abs() < 1e-14, "v2 = {}", v[2]);
    }

    /// Endpoint interpolation property, hand-derived: at u = -1 the value
    /// basis reduces to [1, 0, 0] (only H00 survives) and at u = +1 to
    /// [0, 1, 0] (only H01 survives; W(+-1) = 0 kills the Jacobi part).
    #[test]
    fn hermit_jacobi_endpoint_interpolation_hand() {
        let base = HermitJacobi::new(2, GeomAbsShape::C0);

        let mut vm1 = [0.0f64; 3];
        base.d0(-1.0, &mut vm1);
        assert!((vm1[0] - 1.0).abs() < 1e-14);
        assert!(vm1[1].abs() < 1e-14);
        assert!(vm1[2].abs() < 1e-14);

        let mut vp1 = [0.0f64; 3];
        base.d0(1.0, &mut vp1);
        assert!(vp1[0].abs() < 1e-14);
        assert!((vp1[1] - 1.0).abs() < 1e-14);
        assert!(vp1[2].abs() < 1e-14);
    }
}
