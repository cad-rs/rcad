//! OCCT math_Gauss (TKMath) — Gaussian elimination LU solver.
//!
//! 1:1 translation of `math_Gauss.cxx` + the `LU_Decompose` / `LU_Solve`
//! routines of `math_Recipes.cxx`.

use super::{MatD, VecD};

/// OCCT math_Status_SingularMatrix.
const STATUS_SINGULAR_MATRIX: i32 = 1;
/// OCCT TINY value used by math_Gauss default construction.
const TINY: f64 = 1.0e-20;

/// OCCT math_Gauss: LU factorization of a square matrix with partial pivoting.
#[derive(Debug, Clone)]
pub struct MathGauss {
    lu: MatD,
    index: Vec<i32>,
    d: f64,
    done: bool,
}

impl MathGauss {
    /// OCCT math_Gauss(A, MinPivot) — factorizes A in place into LU.
    pub fn new(a: &MatD) -> Self {
        Self::with_min_pivot(a, TINY)
    }

    /// OCCT math_Gauss(A, MinPivot) with an explicit pivot threshold.
    pub fn with_min_pivot(a: &MatD, min_pivot: f64) -> Self {
        let n = a.n_rows();
        let mut lu = MatD::new(n, a.n_cols());
        for i in 1..=n {
            for j in 1..=a.n_cols() {
                lu.set(i, j, a.get(i, j));
            }
        }
        let mut index = vec![0i32; n];
        let mut d = 0.0f64;
        let error = lu_decompose(&mut lu, &mut index, &mut d, min_pivot);
        MathGauss {
            lu,
            index,
            d,
            done: error == 0,
        }
    }

    /// OCCT math_Gauss::IsDone.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT math_Gauss::Solve(X) — solves in place (X holds B on entry).
    pub fn solve(&self, x: &mut VecD) {
        assert!(self.done, "math_Gauss::Solve - not done");
        lu_solve(&self.lu, &self.index, x);
    }

    /// OCCT math_Gauss::Determinant.
    /// OCCT math_Gauss::Invert(Inv) (math_Gauss.cxx L71-98) — fills a
    /// normalized matrix with the inverse of the LU-decomposed matrix
    /// (column-by-column LU_Solve on unit vectors).
    pub fn invert(&self) -> MatD {
        assert!(self.done, "StdFail_NotDone: math_Gauss::Invert");
        let n = self.lu.n_rows();
        let mut inv = MatD::new(n, self.lu.n_cols());
        for j in 1..=n {
            let mut column = VecD::new(n);
            for i in 1..=n {
                column.set(i, 0.0);
            }
            column.set(j, 1.0);
            lu_solve(&self.lu, &self.index, &mut column);
            for i in 1..=self.lu.n_rows() {
                inv.set(i, j, column.get(i));
            }
        }
        inv
    }

    pub fn determinant(&self) -> f64 {
        assert!(self.done, "math_Gauss::Determinant - not done");
        let mut result = self.d;
        for j in 1..=self.lu.n_rows() {
            result *= self.lu.get(j, j);
        }
        result
    }
}

/// OCCT math_Recipes LU_Decompose (math_Recipes.cxx L298-...): Crout-Doolittle
/// LU with implicit partial pivoting; returns 0 on success.
fn lu_decompose(a: &mut MatD, indx: &mut [i32], d: &mut f64, tiny: f64) -> i32 {
    let n = a.n_rows();
    let mut vv = VecD::new(n);
    *d = 1.0;

    for i in 1..=n {
        let mut big = 0.0f64;
        for j in 1..=n {
            let temp = a.get(i, j).abs();
            if temp > big {
                big = temp;
            }
        }
        if big <= tiny {
            return STATUS_SINGULAR_MATRIX;
        }
        vv.set(i, 1.0 / big);
    }

    for j in 1..=n {
        let mut imax = 0usize;
        for i in 1..j {
            let mut sum = a.get(i, j);
            for k in 1..i {
                sum -= a.get(i, k) * a.get(k, j);
            }
            a.set(i, j, sum);
        }
        let mut big = 0.0f64;
        for i in j..=n {
            let mut sum = a.get(i, j);
            for k in 1..j {
                sum -= a.get(i, k) * a.get(k, j);
            }
            a.set(i, j, sum);
            // Comparison made so imax updates even for NaN/Inf (OCCT #25559).
            let dum = vv.get(i) * sum.abs();
            if dum < big {
                continue;
            }
            big = dum;
            imax = i;
        }
        if j as usize != imax {
            for k in 1..=n {
                let dum = a.get(imax, k);
                a.set(imax, k, a.get(j, k));
                a.set(j, k, dum);
            }
            *d = -*d;
            vv.set(j as usize, vv.get(imax));
        }
        indx[j - 1] = imax as i32;
        if a.get(j, j).abs() <= tiny {
            return STATUS_SINGULAR_MATRIX;
        }
        if j != n {
            let dum = 1.0 / a.get(j, j);
            for i in (j + 1)..=n {
                a.set(i, j, a.get(i, j) * dum);
            }
        }
    }

    0
}

/// OCCT math_Recipes LU_Solve — forward/back substitution with pivot indexing.
fn lu_solve(a: &MatD, indx: &[i32], b: &mut VecD) {
    let n = a.n_rows();
    // OCCT: nblow = b.Lower() - 1 == 0 for 1-based VecD.
    let mut ii = 0usize;
    for i in 1..=n {
        let ip = indx[i - 1] as usize;
        let mut sum = b.get(ip);
        b.set(ip, b.get(i));
        if ii != 0 {
            for j in ii..i {
                sum -= a.get(i, j) * b.get(j);
            }
        } else if sum != 0.0 {
            ii = i;
        }
        b.set(i, sum);
    }
    for i in (1..=n).rev() {
        let mut sum = b.get(i);
        for j in (i + 1)..=n {
            sum -= a.get(i, j) * b.get(j);
        }
        b.set(i, sum / a.get(i, i));
    }
}
