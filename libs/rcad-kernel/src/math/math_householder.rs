// OCCT math_Householder (TKMath/math/math_Householder.hxx/.cxx/.lxx) — 1:1
// Rust translation: least-square resolution of A.X = B via Householder QR
// reflections (the solver of AppParCurves_LeastSquare.gxx Perform's
// NoConstraint branch).
//
// rcad note: rcad's MatD is a 1-based math_Matrix, so LowerRow/LowerCol of A
// are always 1; the submatrix-window constructor still ports verbatim since
// Perform indexes A through mylowerArow/mylowerAcol offsets.

use super::{MatD, VecD};

/// OCCT math_Householder.
#[derive(Debug, Clone)]
pub struct Householder {
    /// OCCT Sol.
    sol: MatD,
    /// OCCT Q.
    q: MatD,
    /// OCCT Done.
    done: bool,
    /// OCCT mylowerArow / myupperArow / mylowerAcol / myupperAcol.
    mylowerarow: usize,
    _myupperarow: usize,
    myloweracol: usize,
    _myupperacol: usize,
}

impl Householder {
    /// OCCT math_Householder(A, B, EPS = 1.0e-20) — the matrix-B constructor
    /// (cxx L28-40).
    pub fn new(a: &MatD, b: &MatD, eps: f64) -> Self {
        let mut h = Householder {
            sol: MatD::new(a.n_cols(), 1),
            q: MatD::new(a.n_rows(), a.n_cols()),
            done: false,
            mylowerarow: 1,
            _myupperarow: a.n_rows(),
            myloweracol: 1,
            _myupperacol: a.n_cols(),
        };
        h.perform(a, b, eps);
        h
    }

    /// OCCT math_Householder(A, B, lowerArow, upperArow, lowerAcol, upperAcol,
    /// EPS = 1.0e-20) — the submatrix-window constructor (cxx L42-58).
    #[allow(clippy::too_many_arguments)]
    pub fn with_window(
        a: &MatD,
        b: &MatD,
        lowerarow: usize,
        upperarow: usize,
        loweracol: usize,
        upperacol: usize,
        eps: f64,
    ) -> Self {
        let mut h = Householder {
            sol: MatD::new(upperacol - loweracol + 1, b.n_cols()),
            q: MatD::new(upperarow - lowerarow + 1, upperacol - loweracol + 1),
            done: false,
            mylowerarow: lowerarow,
            _myupperarow: upperarow,
            myloweracol: loweracol,
            _myupperacol: upperacol,
        };
        h.perform(a, b, eps);
        h
    }

    /// OCCT math_Householder(A, B, EPS = 1.0e-20) — the vector-B constructor
    /// (cxx L16-26): B1(1, B.Length(), 1, 1) with SetCol(1, B).
    pub fn new_vector_b(a: &MatD, b: &VecD, eps: f64) -> Self {
        let mut b1 = MatD::new(b.len(), 1);
        b1.set_col(1, b);
        Householder::new(a, &b1, eps)
    }

    /// OCCT Perform(A, B, EPS) (cxx L60-186).
    fn perform(&mut self, a: &MatD, b: &MatD, eps: f64) {
        let n = self.q.n_cols(); // OCCT n = Q.ColNumber().
        let l = self.q.n_rows(); // OCCT l = Q.RowNumber().
        let m = b.n_cols(); // OCCT m = B.ColNumber().
        let mut b2 = MatD::new(l, m);
        let lbrow = 1; // OCCT lbrow = B.LowerRow().
        for i in 1..=l {
            for j in 1..=n {
                let v = a.get(i + self.mylowerarow - 1, j + self.myloweracol - 1);
                self.q.set(i, j, v);
            }
            for j in 1..=m {
                let v = b.get(i + lbrow - 1, j);
                b2.set(i, j, v);
            }
        }
        // Standard_DimensionError_Raise_if(l != B.RowNumber() || n > l, " ").
        assert!(
            l == b.n_rows() && n <= l,
            "Standard_DimensionError: math_Householder::Perform"
        );

        // Traitement de chaque colonne de A:
        let mut i = 1usize;
        while i <= n {
            let mut h = 0.0;
            for k in i..=l {
                let qki = self.q.get(k, i);
                h += qki * qki; // = ||a||*||a||     = EUAI
            }
            let f = self.q.get(i, i); // = a1              = AII
            let g = if f < 1.0e-15 { h.sqrt() } else { -h.sqrt() };
            if g.abs() <= eps {
                self.done = false;
                return;
            }
            h -= f * g; // = (v*v)/2         = C1
            let alfaii = g - f; // = v               = ALFAII
            for j in (i + 1)..=n {
                let mut scale = 0.0;
                for k in i..=l {
                    scale += self.q.get(k, i) * self.q.get(k, j); // = SCAL
                }
                let cj = (g * self.q.get(i, j) - scale) / h;
                let v = self.q.get(i, j) - alfaii * cj;
                self.q.set(i, j, v);
                for k in (i + 1)..=l {
                    let w = self.q.get(k, j) + cj * self.q.get(k, i);
                    self.q.set(k, j, w);
                }
            }
            // Modification de B:
            for j in 1..=m {
                let mut scale = self.q.get(i, i) * b2.get(i, j);
                for k in (i + 1)..=l {
                    scale += self.q.get(k, i) * b2.get(k, j);
                }
                let cj = (g * b2.get(i, j) - scale) / h;
                let v = b2.get(i, j) - cj * alfaii;
                b2.set(i, j, v);
                for k in (i + 1)..=l {
                    let w = b2.get(k, j) + cj * self.q.get(k, i);
                    b2.set(k, j, w);
                }
            }
            self.q.set(i, i, g);
            i += 1;
        }

        // Remontee:
        for j in 1..=m {
            let v = b2.get(n, j) / self.q.get(n, n);
            self.sol.set(n, j, v);
            let mut i = n as i64 - 1;
            while i >= 1 {
                let iu = i as usize;
                let mut scale = 0.0;
                for k in (iu + 1)..=n {
                    scale += self.q.get(iu, k) * self.sol.get(k, j);
                }
                let v = (b2.get(iu, j) - scale) / self.q.get(iu, iu);
                self.sol.set(iu, j, v);
                i -= 1;
            }
        }
        self.done = true;
    }

    /// OCCT IsDone() (lxx).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Value(sol, Index = 1) (lxx) — the least square solution for the
    /// column Index of B.
    pub fn value(&self, sol: &mut VecD, index: usize) {
        assert!(self.done, "StdFail_NotDone: math_Householder::Value");
        assert!(
            index >= 1 && index <= self.sol.n_cols(),
            "Standard_OutOfRange: math_Householder::Value"
        );
        *sol = self.sol.col(index);
    }

    /// OCCT AllValues() (lxx) — the matrix of all the solutions of A.X = B.
    pub fn all_values(&self) -> &MatD {
        assert!(self.done, "StdFail_NotDone: math_Householder::AllValues");
        &self.sol
    }

    /// OCCT Dump(o) (cxx L188-198).
    pub fn dump(&self) {
        print!("math_Householder ");
        if self.done {
            println!(" Status = Done ");
        } else {
            println!("Status = not Done ");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // TEST(MathHouseholderTest, ExactlyDeterminedSystem) — OCCT
    // math_Householder_Test.cxx: aA(1..3,1..3) * x = b must verify A*x = b.
    #[test]
    fn exactly_determined_system() {
        let mut a = MatD::new(3, 3);
        a.set(1, 1, 1.0);
        a.set(1, 2, 2.0);
        a.set(1, 3, 3.0);
        a.set(2, 1, 4.0);
        a.set(2, 2, 5.0);
        a.set(2, 3, 6.0);
        a.set(3, 1, 7.0);
        a.set(3, 2, 8.0);
        a.set(3, 3, 10.0); // Note: 9 would make it singular
        let mut b = VecD::new(3);
        b.set(1, 14.0);
        b.set(2, 32.0);
        b.set(3, 55.0); // Should give solution approximately [1, 2, 3]

        let householder = Householder::new_vector_b(&a, &b, 1.0e-20);
        assert!(
            householder.is_done(),
            "Householder should succeed for well-conditioned system"
        );

        let mut sol = VecD::new(3);
        householder.value(&mut sol, 1);
        for i in 1..=3 {
            let mut verify = 0.0;
            for j in 1..=3 {
                verify += a.get(i, j) * sol.get(j);
            }
            assert!(
                (verify - b.get(i)).abs() < 1.0e-10,
                "Solution verification A*x=b ({})",
                i
            );
        }
    }

    // TEST(MathHouseholderTest, OverdeterminedSystem) — least squares on 4x2.
    #[test]
    fn overdetermined_system() {
        let mut a = MatD::new(4, 2);
        a.set(1, 1, 1.0);
        a.set(1, 2, 1.0);
        a.set(2, 1, 1.0);
        a.set(2, 2, 2.0);
        a.set(3, 1, 2.0);
        a.set(3, 2, 1.0);
        a.set(4, 1, 1.0);
        a.set(4, 2, 3.0);
        let mut b = VecD::new(4);
        b.set(1, 2.0);
        b.set(2, 3.1);
        b.set(3, 2.9);
        b.set(4, 4.2);

        let householder = Householder::new_vector_b(&a, &b, 1.0e-20);
        assert!(householder.is_done(), "Householder should handle overdetermined system");
        let sol = householder.all_values();
        // Residuals must be finite and small-ish (least squares fit).
        let mut residual = 0.0;
        for i in 1..=4 {
            let mut v = 0.0;
            for j in 1..=2 {
                v += a.get(i, j) * sol.get(j, 1);
            }
            residual += (v - b.get(i)) * (v - b.get(i));
        }
        assert!(residual < 1.0e-4, "residual should be small, got {}", residual);
    }
}
