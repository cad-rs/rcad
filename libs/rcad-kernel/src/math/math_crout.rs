// OCCT math_Crout (TKMath/math/math_Crout.hxx/.cxx) — 1:1 Rust translation:
// Crout-Doolittle symmetric LU decomposition with the inverse computed
// directly from the L·D·L^T factors (the solver used by math_Uzawa's direct
// resolution branch).

use super::math_matrix::{Matrix, Vector};

/// OCCT math_Crout.
#[derive(Debug, Clone)]
pub struct Crout {
    inv_a: Matrix,
    det: f64,
    done: bool,
}

impl Crout {
    /// OCCT math_Crout(A, MinPivot = 1.0e-20) (math_Crout.cxx L17-118).
    pub fn new(a: &Matrix, min_pivot: f64) -> Self {
        let nctl = a.row_number();
        let lowr = a.lower_row();
        let lowc = a.lower_col();
        // math_NotSquare_Raise_if(Nctl != A.ColNumber(), " ").
        assert!(nctl == a.col_number(), "math_NotSquare: math_Crout");

        let mut l = Matrix::new(1, nctl, 1, nctl);
        let mut diag = Vector::new(1, nctl);
        let mut inv_a = Matrix::new(1, nctl, 1, a.col_number());
        let mut det = 1.0;

        for i in 1..=nctl {
            for j in 1..=(i - 1) {
                let mut scale = 0.0;
                for k in 1..=(j - 1) {
                    scale += l.get(i, k) * l.get(j, k) * diag.get(k);
                }
                let v = (a.get(i + lowr - 1, j + lowc - 1) - scale) / diag.get(j);
                l.set(i, j, v);
            }
            let mut scale = 0.0;
            for k in 1..=(i - 1) {
                scale += l.get(i, k) * l.get(i, k) * diag.get(k);
            }
            let d = a.get(i + lowr - 1, i + lowc - 1) - scale;
            diag.set(i, d);
            det *= d;
            if d.abs() <= min_pivot {
                return Crout {
                    inv_a,
                    det,
                    done: false,
                };
            }
            l.set(i, i, 1.0);
        }
        // Calcul de l inverse de L:
        //==========================
        {
            let v = 1.0 / l.get(1, 1);
            l.set(1, 1, v);
        }
        for i in 2..=nctl {
            for k in 1..=(i - 1) {
                let mut scale = 0.0;
                for j in k..=(i - 1) {
                    scale += l.get(i, j) * l.get(j, k);
                }
                let v = -scale / l.get(i, i);
                l.set(i, k, v);
            }
            let v = 1.0 / l.get(i, i);
            l.set(i, i, v);
        }
        // Calcul de l inverse de Mat:
        //============================
        for j in 1..=nctl {
            let mut scale = l.get(j, j) * l.get(j, j) / diag.get(j);
            for k in (j + 1)..=nctl {
                scale += l.get(k, j) * l.get(k, j) / diag.get(k);
            }
            inv_a.set(j, j, scale);
            for i in (j + 1)..=nctl {
                let mut scale = l.get(i, j) * l.get(i, i) / diag.get(i);
                for k in (i + 1)..=nctl {
                    scale += l.get(k, j) * l.get(k, i) / diag.get(k);
                }
                inv_a.set(i, j, scale);
            }
        }
        Crout {
            inv_a,
            det,
            done: true,
        }
    }

    /// OCCT IsDone() (lxx).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Inverse() (hxx) — the inverse of the input matrix.
    pub fn inverse(&self) -> &Matrix {
        &self.inv_a
    }

    /// OCCT Determinant() — the determinant of the input matrix.
    pub fn determinant(&self) -> f64 {
        self.det
    }

    /// OCCT Solve(B, X) (math_Crout.cxx L120-146) — X = InvA * B.
    pub fn solve(&self, b: &Vector, x: &mut Vector) {
        assert!(self.done, "StdFail_NotDone: math_Crout::Solve");
        let n = self.inv_a.row_number();
        // Standard_DimensionError_Raise_if((B.Length() != InvA.RowNumber())
        //     || (X.Length() != B.Length()), " ").
        assert!(
            b.length() == n && x.length() == b.length(),
            "Standard_DimensionError: math_Crout::Solve"
        );
        // OCCT lowb = B.Lower(), lowx = X.Lower() — handled by Vector bounds.
        for i in 1..=n {
            let mut v = self.inv_a.get(i, 1) * b.get(1);
            for j in 2..=i {
                v += self.inv_a.get(i, j) * b.get(j);
            }
            for j in (i + 1)..=n {
                v += self.inv_a.get(j, i) * b.get(j);
            }
            x.set(i, v);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // math_Crout_Test.cxx (OCCT GTests): invert a 3x3 system and verify
    // A * InvA = I, plus Solve against the known solution.
    #[test]
    fn inverse_and_solve_3x3() {
        let mut a = Matrix::new(1, 3, 1, 3);
        a.set(1, 1, 4.0);
        a.set(1, 2, 1.0);
        a.set(1, 3, 0.0);
        a.set(2, 1, 1.0);
        a.set(2, 2, 4.0);
        a.set(2, 3, 1.0);
        a.set(3, 1, 0.0);
        a.set(3, 2, 1.0);
        a.set(3, 3, 4.0);

        let crout = Crout::new(&a, 1.0e-20);
        assert!(crout.is_done());

        let inv = crout.inverse();
        // A * InvA = I. InvA is filled in its lower triangle only (the
        // matrix is symmetric; OCCT's Solve/Uzawa consume it through
        // symmetric access), so read it symmetrically here too.
        for i in 1..=3 {
            for j in 1..=3 {
                let mut v = 0.0;
                for k in 1..=3 {
                    let ij = if k >= j { inv.get(k, j) } else { inv.get(j, k) };
                    v += a.get(i, k) * ij;
                }
                let expect = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (v - expect).abs() < 1.0e-10,
                    "A*InvA({},{}) = {}",
                    i,
                    j,
                    v
                );
            }
        }

        // Solve A x = b with b = (5, 6, 5) -> x = (1, 1, 1).
        let b = Vector::new_init(1, 3, 0.0);
        let mut b = b;
        b.set(1, 5.0);
        b.set(2, 6.0);
        b.set(3, 5.0);
        let mut x = Vector::new(1, 3);
        crout.solve(&b, &mut x);
        assert!((x.get(1) - 1.0).abs() < 1.0e-10);
        assert!((x.get(2) - 1.0).abs() < 1.0e-10);
        assert!((x.get(3) - 1.0).abs() < 1.0e-10);
    }
}
