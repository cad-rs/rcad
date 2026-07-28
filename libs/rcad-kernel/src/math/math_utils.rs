//! OCCT math_utils — legacy forwarding module.
//!
//! All functions have been moved to OCCT Math*-package sub-modules:
//! - root/    — MathRoot (newton_raphson, bisection, secant, biss_newton, trig_roots, etc.)
//! - poly/    — MathPoly (solve_linear/quadratic/cubic/quartic, laguerre_roots, poly_eval)
//! - opt/     — MathOpt (bfgs, frpr, newton_minimize, powell, brent, golden_section, etc.)
//! - lin/     — MathLin (svd, gauss, crout, householder, eigenvalues, matrix, least_squares)
//! - integ/   — MathInteg (simpson, gaussian_quadrature)
//! - sys/     — MathSys (newton_2d, newton_3d)
//! - gprop/plate/ — thin_plate_spline, evaluate_tps, build_plate_surface
//!
//! This file is kept for backward compatibility during migration.
//! Downstream code should import from the new sub-modules directly.

pub use crate::math::root::{
    newton_raphson, bisection, secant, biss_newton,
    trig_roots, trig_roots_sin_only, trig_roots_cos_sin,
    find_roots_in, bracket_root,
};
pub use crate::math::math_poly::{
    solve_linear, solve_quadratic, solve_cubic, solve_quartic,
    laguerre_roots, poly_eval,
};
pub use crate::math::opt::{
    bfgs_minimize, frpr_minimize, newton_minimize, powell_minimize,
    golden_section_min, golden_section_max, brent_minimize,
    glob_opt_min, pso_minimize, lm_solve,
};
pub use crate::math::lin::{
    svd_solve_3x3, eigenvalues_2x2, eigenvalues_3x3,
    inverse_3x3, determinant_3x3,
    solve_linear_system, crout_solve, householder_solve,
    least_squares_linear,
};
pub use crate::math::integ::{simpson_integrate, gaussian_quadrature};
pub use crate::math::sys::{newton_2d, newton_3d};
pub use crate::math::gprop::plate::{
    thin_plate_spline, evaluate_tps, build_plate_surface,
};
