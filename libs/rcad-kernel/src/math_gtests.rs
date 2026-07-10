//! OCCT-aligned math_* GTest translations.
//!
//! OCCT source: src/FoundationClasses/TKMath/GTests/
//!
//! Files translated:
//!   math_SVD_Test.cxx                 — SVD decomposition and solve
//!   math_GaussLeastSquare_Test.cxx    — Least squares fitting
//!   math_BFGS_Test.cxx                — BFGS optimization
//!   math_FRPR_Test.cxx                — Fletcher-Reeves Polak-Ribiere
//!   math_NewtonMinimum_Test.cxx       — Newton optimization with Hessian
//!   math_FunctionRoot_Test.cxx        — Newton-Raphson, bisection, secant
//!   math_DirectPolynomialRoots_Test   — solve_quadratic, solve_cubic, solve_quartic
//!   math_NewtonFunctionSetRoot_Test   — newton_2d, newton_3d
//!   math_Integration_Test             — simpson_integrate, gaussian_quadrature
//!   math_Matrix_Test                  — determinant, eigenvalues, inverse
//!   math_Vector_Test                  — vector norm, dot, cross
//!   math_BracketMinimum_Test          — golden_section_min (bracket-then-scan)
//!   MathOpt_1D_Test                   — golden_section_min/max
//!   MathPoly_Test                     — polynomial evaluation
//!   MathLin_EigenSearch_Test          — eigenvalue search
//!   math_Gauss_Test                   — Gaussian elimination

use glam::{DMat2, DMat3, DVec2, DVec3};

const TOL: f64 = 1e-10;

// =============================================================================
// math_FunctionRoot_Test.cxx — Newton-Raphson, bisection, secant
// =============================================================================

#[cfg(test)]
mod function_root_tests {
    use super::*;

    fn f_quadratic(x: f64) -> f64 { x * x - 4.0 }          // roots at ±2
    fn df_quadratic(x: f64) -> f64 { 2.0 * x }
    fn f_sin(x: f64) -> f64 { x.sin() - 0.5 }               // root at π/6
    fn d_sin(x: f64) -> f64 { x.cos() }

    #[test]
    fn newton_raphson_quadratic() {
        let root = crate::math_utils::newton_raphson(f_quadratic, df_quadratic, 3.0, 1e-12, 100);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0).abs() < 1e-6);
    }

    #[test]
    fn newton_raphson_sin() {
        let root = crate::math_utils::newton_raphson(f_sin, d_sin, 1.0, 1e-12, 100);
        assert!(root.is_some());
        assert!((root.unwrap() - std::f64::consts::FRAC_PI_6).abs() < 1e-6);
    }

    #[test]
    fn bisection_quadratic() {
        let root = crate::math_utils::bisection(f_quadratic, 0.0, 5.0, 1e-12);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0).abs() < 1e-6);
    }

    #[test]
    fn bisection_no_root_in_interval() {
        // f(x) = x^2 - 4 has no root in [0, 1]
        let root = crate::math_utils::bisection(f_quadratic, 0.0, 1.0, 1e-12);
        assert!(root.is_none());
    }

    #[test]
    fn secant_quadratic() {
        let root = crate::math_utils::secant(f_quadratic, 1.0, 3.0, 1e-12);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0).abs() < 1e-6);
    }

    #[test]
    fn newton_raphson_no_convergence() {
        // f(x) = x^2 + 1 has no real root
        fn f_no_root(x: f64) -> f64 { x * x + 1.0 }
        fn df_no_root(x: f64) -> f64 { 2.0 * x }
        let root = crate::math_utils::newton_raphson(f_no_root, df_no_root, 0.0, 1e-12, 10);
        assert!(root.is_none());
    }
}

// =============================================================================
// math_NewtonFunctionSetRoot_Test.cxx — newton_2d, newton_3d
// =============================================================================

#[cfg(test)]
mod newton_set_root_tests {
    use super::*;

    #[test]
    fn newton_2d_simple_intersection() {
        // Two lines: x + y = 3, x - y = 1 → solution: (2, 1)
        fn f(v: DVec2) -> DVec2 {
            DVec2::new(v.x + v.y - 3.0, v.x - v.y - 1.0)
        }
        fn jacobian(_v: DVec2) -> DMat2 {
            DMat2::from_cols(DVec2::new(1.0, 1.0), DVec2::new(1.0, -1.0))
        }
        let sol = crate::math_utils::newton_2d(f, jacobian, DVec2::new(0.0, 0.0), 1e-10);
        assert!(sol.is_some());
        let s = sol.unwrap();
        assert!((s.x - 2.0).abs() < 1e-6);
        assert!((s.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn newton_3d_linear_system() {
        // 3x3 linear system: Ax = b where A is identity, b = (1,2,3)
        fn f(v: DVec3) -> DVec3 {
            DVec3::new(v.x - 1.0, v.y - 2.0, v.z - 3.0)
        }
        fn jacobian(_v: DVec3) -> DMat3 {
            DMat3::IDENTITY
        }
        let sol = crate::math_utils::newton_3d(f, jacobian, DVec3::ZERO, 1e-10);
        assert!(sol.is_some());
        let s = sol.unwrap();
        assert!((s - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-6);
    }
}

// =============================================================================
// math_DirectPolynomialRoots_Test.cxx — Polynomial solvers
// =============================================================================

#[cfg(test)]
mod polynomial_roots_tests {
    use super::*;

    #[test]
    fn solve_linear_basic() {
        let root = crate::math_utils::solve_linear(2.0, -4.0);
        assert!((root.unwrap() - 2.0).abs() < TOL);
    }

    #[test]
    fn solve_linear_no_solution() {
        let root = crate::math_utils::solve_linear(0.0, 5.0);
        assert!(root.is_none());
    }

    #[test]
    fn solve_quadratic_two_roots() {
        let roots = crate::math_utils::solve_quadratic(1.0, 0.0, -4.0); // x^2 - 4 = 0
        assert_eq!(roots.len(), 2);
        assert!((roots[0] + 2.0).abs() < TOL); // sorted: -2, 2
        assert!((roots[1] - 2.0).abs() < TOL);
    }

    #[test]
    fn solve_quadratic_no_real_roots() {
        let roots = crate::math_utils::solve_quadratic(1.0, 0.0, 1.0); // x^2 + 1 = 0
        assert!(roots.is_empty());
    }

    #[test]
    fn solve_quadratic_double_root() {
        let roots = crate::math_utils::solve_quadratic(1.0, -4.0, 4.0); // (x-2)^2 = 0
        assert_eq!(roots.len(), 1);
        assert!((roots[0] - 2.0).abs() < TOL);
    }

    #[test]
    fn solve_cubic_three_roots() {
        let roots = crate::math_utils::solve_cubic(1.0, -6.0, 11.0, -6.0); // (x-1)(x-2)(x-3)
        assert_eq!(roots.len(), 3);
        assert!((roots[0] - 1.0).abs() < TOL);
        assert!((roots[1] - 2.0).abs() < TOL);
        assert!((roots[2] - 3.0).abs() < TOL);
    }

    #[test]
    fn solve_cubic_single_root() {
        // x^3 - 1 = 0 has one real root (x=1)
        let roots = crate::math_utils::solve_cubic(1.0, 0.0, 0.0, -1.0);
        assert!(roots.len() >= 1);
        assert!((roots[0] - 1.0).abs() < TOL);
    }

    #[test]
    fn solve_quartic_four_roots() {
        // x^4 - 10x^2 + 9 = 0 → roots: ±1, ±3
        let roots = crate::math_utils::solve_quartic(1.0, 0.0, -10.0, 0.0, 9.0);
        assert_eq!(roots.len(), 4);
        assert!((roots[0] + 3.0).abs() < TOL);
        assert!((roots[1] + 1.0).abs() < TOL);
        assert!((roots[2] - 1.0).abs() < TOL);
        assert!((roots[3] - 3.0).abs() < TOL);
    }
}

// =============================================================================
// math_Integration_Test.cxx + MathInteg_Test.cxx — Numerical integration
// =============================================================================

#[cfg(test)]
mod integration_tests {
    use super::*;

    fn f_x_squared(x: f64) -> f64 { x * x }
    fn f_sin(x: f64) -> f64 { x.sin() }

    #[test]
    fn simpson_x_squared() {
        // ∫₀¹ x² dx = 1/3
        let result = crate::math_utils::simpson_integrate(f_x_squared, 0.0, 1.0, 100);
        assert!((result - 1.0 / 3.0).abs() < 1e-4);
    }

    #[test]
    fn simpson_sin() {
        // ∫₀^π sin(x) dx = 2
        let result = crate::math_utils::simpson_integrate(f_sin, 0.0, std::f64::consts::PI, 100);
        assert!((result - 2.0).abs() < 1e-4);
    }

    #[test]
    fn gaussian_quadrature_x_squared() {
        let result = crate::math_utils::gaussian_quadrature(f_x_squared, 0.0, 1.0, 5);
        assert!((result - 1.0 / 3.0).abs() < 1e-4);
    }

    #[test]
    fn gaussian_quadrature_sin() {
        let result = crate::math_utils::gaussian_quadrature(f_sin, 0.0, std::f64::consts::PI, 5);
        assert!((result - 2.0).abs() < 1e-4);
    }
}

// =============================================================================
// math_Matrix_Test.cxx + MathLin_EigenSearch_Test.cxx — Matrix operations
// =============================================================================

#[cfg(test)]
mod matrix_tests {
    use super::*;

    #[test]
    fn determinant_identity() {
        assert!((crate::math_utils::determinant_3x3(DMat3::IDENTITY) - 1.0).abs() < TOL);
    }

    #[test]
    fn determinant_scale() {
        let m = DMat3::from_diagonal(DVec3::new(2.0, 3.0, 4.0));
        assert!((crate::math_utils::determinant_3x3(m) - 24.0).abs() < TOL);
    }

    #[test]
    fn inverse_identity() {
        let inv = crate::math_utils::inverse_3x3(DMat3::IDENTITY);
        assert!(inv.is_some());
        assert!((inv.unwrap() - DMat3::IDENTITY).x_axis.length() < TOL);
    }

    #[test]
    fn inverse_scale_then_roundtrip() {
        let m = DMat3::from_diagonal(DVec3::new(2.0, 3.0, 4.0));
        let inv = crate::math_utils::inverse_3x3(m).unwrap();
        let v = DVec3::new(1.0, 2.0, 3.0);
        assert!((inv * (m * v) - v).length() < TOL);
    }

    #[test]
    fn eigenvalues_2x2_diagonal() {
        let m = DMat2::from_diagonal(DVec2::new(2.0, 3.0));
        let (e1, e2) = crate::math_utils::eigenvalues_2x2(m);
        assert!(e1.is_finite() && e2.is_finite());
        assert!((e1 * e2 - m.determinant()).abs() < TOL);
    }

    #[test]
    fn eigenvalues_3x3_diagonal() {
        let m = DMat3::from_diagonal(DVec3::new(1.0, 2.0, 3.0));
        let (e1, e2, e3) = crate::math_utils::eigenvalues_3x3(m);
        assert!(e1.is_finite() && e2.is_finite() && e3.is_finite());
        assert!((e1 * e2 * e3 - m.determinant()).abs() < TOL);
    }
}

// =============================================================================
// math_Vector_Test.cxx — Vector operations (DVec3)
// =============================================================================

#[cfg(test)]
mod math_vector_tests {
    use super::*;

    #[test]
    fn vector_length() {
        assert!((DVec3::new(3.0, 4.0, 0.0).length() - 5.0).abs() < TOL);
    }

    #[test]
    fn vector_normalize() {
        let v = DVec3::new(3.0, 4.0, 0.0).normalize();
        assert!((v.length() - 1.0).abs() < TOL);
    }

    #[test]
    fn vector_dot() {
        assert!((DVec3::new(1.0, 2.0, 3.0).dot(DVec3::new(4.0, 5.0, 6.0)) - 32.0).abs() < TOL);
    }

    #[test]
    fn vector_cross() {
        assert!((DVec3::X.cross(DVec3::Y) - DVec3::Z).length() < TOL);
    }
}

// =============================================================================
// MathOpt_1D_Test.cxx / math_BracketMinimum_Test — Golden section optimization
// =============================================================================

#[cfg(test)]
mod optimization_tests {
    use super::*;

    #[test]
    fn golden_section_min_quadratic() {
        // f(x) = (x-2)^2 + 1, minimum at x=2
        let xmin = crate::math_utils::golden_section_min(|x| (x - 2.0) * (x - 2.0) + 1.0, 0.0, 5.0, 1e-8);
        assert!((xmin - 2.0).abs() < 1e-6);
    }

    #[test]
    fn golden_section_min_sin() {
        // sin(x) minimum near 3π/2 in [2, 5]
        let xmin = crate::math_utils::golden_section_min(|x| x.sin(), 2.0, 5.0, 1e-8);
        assert!((xmin - 3.0 * std::f64::consts::FRAC_PI_2).abs() < 1e-6);
    }

    #[test]
    fn golden_section_max_negative_quadratic() {
        // f(x) = -(x-3)^2 + 5, maximum at x=3
        let xmax = crate::math_utils::golden_section_max(|x| -(x - 3.0) * (x - 3.0) + 5.0, 0.0, 6.0, 1e-8);
        assert!((xmax - 3.0).abs() < 1e-6);
    }
}

// =============================================================================
// MathPoly_Test.cxx — Polynomial evaluation utilities
// =============================================================================

#[cfg(test)]
mod polynomial_eval_tests {
    use super::*;

    /// Horner's method polynomial evaluation: a₀ + a₁x + a₂x² + ... + aₙxⁿ
    fn poly_eval(coeffs: &[f64], x: f64) -> f64 {
        coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
    }

    #[test]
    fn polynomial_constant() {
        assert!((poly_eval(&[5.0], 10.0) - 5.0).abs() < TOL);
    }

    #[test]
    fn polynomial_linear() {
        // 2x + 3 at x = 4 → 11
        assert!((poly_eval(&[3.0, 2.0], 4.0) - 11.0).abs() < TOL);
    }

    #[test]
    fn polynomial_quadratic() {
        // x² + 2x + 1 at x = 3 → 16
        assert!((poly_eval(&[1.0, 2.0, 1.0], 3.0) - 16.0).abs() < TOL);
    }

    #[test]
    fn polynomial_cubic() {
        // x³ - 6x² + 11x - 6 at x = 2 → 0
        let r = poly_eval(&[-6.0, 11.0, -6.0, 1.0], 2.0);
        assert!(r.abs() < TOL);
    }
}

// =============================================================================
// math_Gauss_Test.cxx — Gaussian elimination (via DMat3 inverse)
// =============================================================================

#[cfg(test)]
mod gauss_tests {
    use super::*;

    #[test]
    fn matrix_solve_3x3_identity() {
        // A * x = b, A = I, b = (1,2,3) → x = (1,2,3)
        let a = DMat3::IDENTITY;
        let inv = crate::math_utils::inverse_3x3(a).unwrap();
        let x = inv * DVec3::new(1.0, 2.0, 3.0);
        assert!((x - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn matrix_solve_3x3_diagonal() {
        let a = DMat3::from_diagonal(DVec3::new(2.0, 3.0, 4.0));
        let inv = crate::math_utils::inverse_3x3(a).unwrap();
        let x = inv * DVec3::new(6.0, 12.0, 20.0);
        assert!((x - DVec3::new(3.0, 4.0, 5.0)).length() < TOL);
    }
}

// =============================================================================
// math_SVD_Test.cxx — SVD decomposition and solve
// =============================================================================

#[cfg(test)]
mod svd_tests {
    use super::*;

    #[test]
    fn svd_well_conditioned_3x3() {
        // Matrix: [[2,1,0],[1,2,1],[0,1,2]]
        let a = DMat3::from_cols(
            DVec3::new(2.0, 1.0, 0.0),
            DVec3::new(1.0, 2.0, 1.0),
            DVec3::new(0.0, 1.0, 2.0),
        );
        let b = DVec3::new(6.0, 9.0, 8.0);

        // Solve via SVD
        let x = crate::math_utils::svd_solve_3x3(a, b);
        assert!(x.is_some(), "SVD solve should succeed");

        // Verify A*x = b
        let ax = a * x.unwrap();
        assert!((ax - b).length() < 1e-4, "A*x should equal b");
    }

    #[test]
    fn svd_identity_3x3() {
        let a = DMat3::IDENTITY;
        let b = DVec3::new(1.0, 2.0, 3.0);
        let x = crate::math_utils::svd_solve_3x3(a, b);
        assert!(x.is_some());
        assert!((x.unwrap() - b).length() < TOL);
    }

    #[test]
    fn svd_singular_matrix() {
        // Singular: [[1,1],[1,1]]
        let a = DMat3::from_cols(
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
        );
        let b = DVec3::new(2.0, 2.0, 1.0);
        let x = crate::math_utils::svd_solve_3x3(a, b);
        assert!(x.is_some(), "SVD should handle singular matrices");
        // Verify A*x ≈ b
        let ax = a * x.unwrap();
        assert!((ax - b).length() < 1e-8);
    }
}

// =============================================================================
// math_GaussLeastSquare_Test.cxx — Least squares fitting
// =============================================================================

#[cfg(test)]
mod least_square_tests {
    use super::*;

    #[test]
    fn least_squares_line_fit() {
        // Points: (1,1), (2,2), (3,3), (4,4) — perfect line y=x
        let x_pts = vec![1.0, 2.0, 3.0, 4.0];
        let y_pts = vec![1.0, 2.0, 3.0, 4.0];
        // Fit y = a + b*x
        let result = crate::math_utils::least_squares_linear(&x_pts, &y_pts);
        assert!(result.is_some(), "LS fit should succeed");
        let (a, b) = result.unwrap();
        assert!((a - 0.0).abs() < 1e-8, "intercept should be 0");
        assert!((b - 1.0).abs() < 1e-8, "slope should be 1");
    }

    #[test]
    fn least_squares_noisy() {
        let x_pts = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y_pts = vec![0.1, 2.1, 3.9, 6.0, 8.1]; // y ≈ 2x
        let result = crate::math_utils::least_squares_linear(&x_pts, &y_pts);
        assert!(result.is_some());
        let (_a, b) = result.unwrap();
        assert!((b - 2.0).abs() < 0.2, "slope should be ~2, got {b}");
    }
}

// =============================================================================
// math_BFGS_Test.cxx / math_FRPR_Test.cxx — Quasi-Newton optimization
// =============================================================================

#[cfg(test)]
mod bfgs_tests {
    use super::*;

    /// Quadratic bowl: f(x,y) = (x-1)^2 + (y-2)^2, minimum at (1,2)
    fn quadratic_bowl_grad(x: &[f64], g: &mut [f64]) -> f64 {
        g[0] = 2.0 * (x[0] - 1.0);
        g[1] = 2.0 * (x[1] - 2.0);
        (x[0] - 1.0).powi(2) + (x[1] - 2.0).powi(2)
    }

    #[test]
    fn bfgs_quadratic_bowl() {
        let x0 = vec![0.0, 0.0];
        let result = crate::math_utils::bfgs_minimize(&x0, quadratic_bowl_grad, 1e-10, 100);
        assert!(result.is_some(), "BFGS should converge on quadratic bowl");
        let x = result.unwrap();
        assert!((x[0] - 1.0).abs() < 1e-6, "x should be 1, got {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-6, "y should be 2, got {}", x[1]);
    }

    #[test]
    fn bfgs_rosenbrock() {
        // Rosenbrock: f(x,y) = 100*(y-x^2)^2 + (1-x)^2, minimum at (1,1)
        fn rosenbrock_grad(x: &[f64], g: &mut [f64]) -> f64 {
            g[0] = -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]);
            g[1] = 200.0 * (x[1] - x[0] * x[0]);
            100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2)
        }
        let x0 = vec![-1.0, 1.0];
        let result = crate::math_utils::bfgs_minimize(&x0, rosenbrock_grad, 1e-8, 200);
        assert!(result.is_some(), "BFGS should converge on Rosenbrock");
        let x = result.unwrap();
        assert!((x[0] - 1.0).abs() < 1e-4, "x should be 1, got {}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-4, "y should be 1, got {}", x[1]);
    }
}

// =============================================================================
// math_NewtonMinimum_Test.cxx — Newton optimization with Hessian
// =============================================================================

#[cfg(test)]
mod newton_min_tests {
    use super::*;

    /// Quadratic: f(x,y) = (x-1)^2 + 2*(y-2)^2, minimum at (1,2)
    fn quad_bowl_hessian(x: &[f64], g: &mut [f64], h: &mut [f64]) -> f64 {
        g[0] = 2.0 * (x[0] - 1.0);
        g[1] = 4.0 * (x[1] - 2.0);
        // Hessian: [[2,0],[0,4]]
        h[0] = 2.0; h[1] = 0.0;
        h[2] = 0.0; h[3] = 4.0;
        (x[0] - 1.0).powi(2) + 2.0 * (x[1] - 2.0).powi(2)
    }

    #[test]
    fn newton_min_quadratic_bowl() {
        let x0 = vec![0.0, 0.0];
        let result = crate::math_utils::newton_minimize(&x0, quad_bowl_hessian, 1e-10, 50);
        assert!(result.is_some(), "Newton should converge on quadratic bowl");
        let x = result.unwrap();
        assert!((x[0] - 1.0).abs() < 1e-8);
        assert!((x[1] - 2.0).abs() < 1e-8);
    }
}
