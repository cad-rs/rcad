//! Newton-Raphson GCS solver with numerical Jacobian.
//!
//! The solver treats the constraint system as F(x) = 0 where x is the vector
//! of free (non-fixed) parameters.  It iterates:
//!
//!   Δx = −J⁺ · F(x)
//!   x  ← x + Δx
//!
//! where J⁺ is the Moore-Penrose pseudo-inverse of the Jacobian, computed via
//! the normal equations (J^T J + λI) Δx = J^T (−F).  The Tikhonov
//! regularisation λ handles under-constrained systems by returning the
//! minimum-norm step.
//!
//! The Jacobian is approximated by central finite differences.

use crate::constraint::Constraint;
use crate::entity::Entity;

/// Result returned by [`crate::Sketch::solve`].
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// `true` if all constraint residuals are below [`RESIDUAL_TOL`].
    pub converged: bool,
    /// RMS residual at termination.
    pub residual: f64,
    /// Number of Newton iterations performed.
    pub iterations: usize,
}

/// Convergence tolerance on the RMS constraint residual.
pub const RESIDUAL_TOL: f64 = 1e-10;
/// Maximum Newton iterations.
const MAX_ITER: usize = 100;
/// Step size for numerical Jacobian (central differences).
const FD_H: f64 = 1e-7;
/// Tikhonov regularisation coefficient.
const LAMBDA: f64 = 1e-10;

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Run the Newton-Raphson solver.
///
/// - `params`: full parameter vector (modified in-place for free params).
/// - `fixed`: mask — `fixed[i] == true` means `params[i]` is held constant.
/// - `entities`: entity metadata (for constraint evaluation).
/// - `constraints`: list of constraints.
pub fn solve(
    params: &mut Vec<f64>,
    fixed: &[bool],
    entities: &[Entity],
    constraints: &[Constraint],
) -> SolveResult {
    // Build index maps: free_params[k] = global param index of the k-th free param.
    let free_params: Vec<usize> = (0..params.len()).filter(|&i| !fixed[i]).collect();
    let n_free = free_params.len();

    // Total number of constraint equations.
    let n_eq: usize = constraints.iter().map(|c| c.equation_count()).sum();

    if n_eq == 0 || n_free == 0 {
        return SolveResult { converged: true, residual: 0.0, iterations: 0 };
    }

    let mut residual;
    let mut iters = 0;

    for _ in 0..MAX_ITER {
        let f = eval_residuals(params, entities, constraints, n_eq);
        residual = rms(&f);
        if residual < RESIDUAL_TOL {
            break;
        }

        // Build numerical Jacobian J (n_eq × n_free)
        let j = numerical_jacobian(params, fixed, &free_params, entities, constraints, n_eq);

        // Solve (J^T J + λI) Δx = J^T (−F)  →  Δx
        let delta = solve_normal_equations(&j, &f, n_free);

        // Update free parameters with step-size damping (halve step if residual grows).
        let mut alpha = 1.0_f64;
        let old_params: Vec<f64> = free_params.iter().map(|&gi| params[gi]).collect();
        for _ in 0..8 {
            for (k, &gi) in free_params.iter().enumerate() {
                params[gi] = old_params[k] + alpha * delta[k];
            }
            let f_new = eval_residuals(params, entities, constraints, n_eq);
            if rms(&f_new) < residual * (1.0 + 1e-4) {
                break;
            }
            alpha *= 0.5;
        }

        iters += 1;
    }

    // Final residual check
    let f_final = eval_residuals(params, entities, constraints, n_eq);
    residual = rms(&f_final);

    SolveResult {
        converged: residual < RESIDUAL_TOL,
        residual,
        iterations: iters,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Residual evaluation
// ─────────────────────────────────────────────────────────────────────────────

pub(crate) fn eval_residuals(
    params: &[f64],
    entities: &[Entity],
    constraints: &[Constraint],
    n_eq: usize,
) -> Vec<f64> {
    let mut f = vec![0.0_f64; n_eq];
    let mut row = 0;
    for c in constraints {
        let cnt = c.equation_count();
        c.residuals(params, entities, &mut f[row..row + cnt]);
        row += cnt;
    }
    f
}

// ─────────────────────────────────────────────────────────────────────────────
// Numerical Jacobian (central differences)
// ─────────────────────────────────────────────────────────────────────────────

fn numerical_jacobian(
    params: &[f64],
    _fixed: &[bool],
    free_params: &[usize],
    entities: &[Entity],
    constraints: &[Constraint],
    n_eq: usize,
) -> Vec<Vec<f64>> {
    let n_free = free_params.len();
    // j[row][col] = ∂F_row / ∂x_col
    let mut j = vec![vec![0.0_f64; n_free]; n_eq];
    let mut p = params.to_vec();

    for (col, &gi) in free_params.iter().enumerate() {
        let orig = p[gi];

        p[gi] = orig + FD_H;
        let f_plus = eval_residuals(&p, entities, constraints, n_eq);

        p[gi] = orig - FD_H;
        let f_minus = eval_residuals(&p, entities, constraints, n_eq);

        p[gi] = orig;

        for row in 0..n_eq {
            j[row][col] = (f_plus[row] - f_minus[row]) / (2.0 * FD_H);
        }
    }
    j
}

// ─────────────────────────────────────────────────────────────────────────────
// Normal equations solver  (J^T J + λI) Δx = J^T (−F)
// ─────────────────────────────────────────────────────────────────────────────

fn solve_normal_equations(j: &[Vec<f64>], f: &[f64], n: usize) -> Vec<f64> {
    let m = j.len();

    // A = J^T J + λI  (n×n)
    let mut a = vec![vec![0.0_f64; n]; n];
    // b = J^T (−F)    (n-vector)
    let mut b = vec![0.0_f64; n];

    for i in 0..m {
        for col in 0..n {
            b[col] -= j[i][col] * f[i];
            for k in 0..n {
                a[col][k] += j[i][col] * j[i][k];
            }
        }
    }
    // Tikhonov regularisation
    for i in 0..n {
        a[i][i] += LAMBDA;
    }

    gaussian_elimination(&mut a, &mut b).unwrap_or_else(|| vec![0.0; n])
}

// ─────────────────────────────────────────────────────────────────────────────
// Gaussian elimination with partial pivoting
// ─────────────────────────────────────────────────────────────────────────────

fn gaussian_elimination(a: &mut Vec<Vec<f64>>, b: &mut Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    debug_assert_eq!(a.len(), n);

    for col in 0..n {
        // Find pivot row (max absolute value in column)
        let pivot_row = (col..n)
            .max_by(|&i, &j| a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap())?;
        a.swap(col, pivot_row);
        b.swap(col, pivot_row);

        let pivot = a[col][col];
        if pivot.abs() < 1e-14 {
            return None;
        }

        // Eliminate below
        for row in (col + 1)..n {
            let factor = a[row][col] / pivot;
            for k in col..n {
                let v = a[col][k] * factor;
                a[row][k] -= v;
            }
            let bv = b[col] * factor;
            b[row] -= bv;
        }
    }

    // Back substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        x[i] = b[i];
        for j in (i + 1)..n {
            let v = a[i][j] * x[j];
            x[i] -= v;
        }
        x[i] /= a[i][i];
    }
    Some(x)
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn rms(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests for the linear algebra primitives
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn debug_fixed_point_inline() {
        use crate::constraint::Constraint;
        use crate::entity::{Entity, EntityKind, PointRef};

        let entities = vec![Entity::new(EntityKind::Point, 0)];
        let constraints = vec![Constraint::Fixed {
            point: PointRef::Point(0),
            x: 0.0,
            y: 0.0,
        }];
        let mut params = vec![3.0_f64, 7.0];
        let fixed = vec![false, false];

        let n_eq = 2usize;
        let f = eval_residuals(&params, &entities, &constraints, n_eq);
        println!("initial f = {:?}", f);

        let free_params = vec![0usize, 1];
        let j = numerical_jacobian(&params, &fixed, &free_params, &entities, &constraints, n_eq);
        println!("J = {:?}", j);

        let delta = solve_normal_equations(&j, &f, 2);
        println!("delta = {:?}", delta);

        let result = solve(&mut params, &fixed, &entities, &constraints);
        println!("params after = {:?}", params);
        println!("converged={} residual={}", result.converged, result.residual);
        assert!(result.converged, "not converged: {}", result.residual);
        assert!((params[0] - 0.0).abs() < 1e-7);
        assert!((params[1] - 0.0).abs() < 1e-7);
    }

    #[test]
    fn gaussian_2x2() {
        // 2x + y = 5, x + 3y = 10  →  x=1, y=3
        let mut a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let mut b = vec![5.0, 10.0];
        let x = gaussian_elimination(&mut a, &mut b).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-10, "x={}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-10, "y={}", x[1]);
    }

    #[test]
    fn gaussian_3x3() {
        // x + y + z = 6, 2x + y = 5, y + 3z = 10  →  x=2, y=1, z=3
        let mut a = vec![
            vec![1.0, 1.0, 1.0],
            vec![2.0, 1.0, 0.0],
            vec![0.0, 1.0, 3.0],
        ];
        let mut b = vec![6.0, 5.0, 10.0];
        let x = gaussian_elimination(&mut a, &mut b).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-9, "x={}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-9, "y={}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-9, "z={}", x[2]);
    }
}
