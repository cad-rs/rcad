// OCCT math_SVD (FoundationClasses/TKMath/math/math_SVD.hxx L33-72,
// math_SVD.cxx L29-116, math_SVD.lxx) — Singular Value Decomposition solver.
//
// SVD implements the solution of a set of N linear equations of M unknowns
// without condition on N or M. For singular or nearly singular matrices SVD
// is a better choice than Gauss or GaussLeastSquare.
//
// The OCCT member fields are (math_SVD.hxx L66-71):
//   bool Done; math_Matrix U; math_Matrix V; math_Vector Diag; int RowA;
// The bound-aware math_Matrix/math_Vector shells are mapped to the kernel
// MatD/VecD forms (architecture difference: MatD/VecD are always 1-based and
// bounds-less, mirroring what math_gauss.rs already does for math_Gauss).

use super::math_recipes::{svd_decompose, svd_solve, MATH_STATUS_OK};
use super::{MatD, VecD};
use crate::core::precision::REAL_FIRST;

/// OCCT math_SVD.
#[derive(Debug, Clone)]
pub struct MathSvd {
    done: bool,
    u: MatD,
    v: MatD,
    diag: VecD,
    row_a: usize,
}

impl MathSvd {
    /// OCCT math_SVD::math_SVD(const math_Matrix& A) (math_SVD.cxx L29-39) —
    /// given an n X m matrix A with n < m, n = m or n > m, performs the
    /// Singular Value Decomposition.
    pub fn new(a: &MatD) -> Self {
        let row_number = a.n_rows(); // OCCT A.RowNumber()
        let col_number = a.n_cols(); // OCCT A.ColNumber()
        // OCCT L30-32: U(1, max(RowNumber, ColNumber), 1, ColNumber),
        //              V(1, ColNumber, 1, ColNumber), Diag(1, ColNumber).
        let mut u = MatD::new(row_number.max(col_number), col_number);
        let mut v = MatD::new(col_number, col_number);
        let mut diag = VecD::new(col_number);
        // OCCT L34: U.Init(0.0);
        for i in 1..=u.n_rows() {
            for j in 1..=u.n_cols() {
                u.set(i, j, 0.0);
            }
        }
        let row_a = row_number; // OCCT L35: RowA = A.RowNumber();
        // OCCT L36: U.Set(1, A.RowNumber(), 1, A.ColNumber(), A);
        for i in 1..=row_number {
            for j in 1..=col_number {
                u.set(i, j, a.get(i, j));
            }
        }
        // OCCT L37-38: int Error = SVD_Decompose(U, Diag, V); Done = Error == 0;
        let error = svd_decompose(&mut u, &mut diag, &mut v);
        MathSvd {
            done: error == MATH_STATUS_OK,
            u,
            v,
            diag,
            row_a,
        }
    }

    /// OCCT math_SVD::IsDone (math_SVD.lxx) — true if the computations are
    /// successful.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT math_VectorBase::Max (math_VectorBase.lxx L162-176) — the 1-based
    /// INDEX of the maximum coefficient (RealFirst seed).
    fn diag_max_index(&self) -> usize {
        let mut i_max = 0usize;
        let mut x = REAL_FIRST;
        for index in 1..=self.diag.len() {
            if self.diag.get(index) > x {
                x = self.diag.get(index);
                i_max = index;
            }
        }
        i_max
    }

    /// OCCT math_SVD::Solve(B, X, Eps) (math_SVD.cxx L41-68) — solves the set
    /// of linear equations A . X = B.
    pub fn solve(&mut self, b: &VecD, x: &mut VecD, eps: f64) {
        // OCCT L43: StdFail_NotDone_Raise_if(!Done, " ");
        assert!(self.done, "StdFail_NotDone: math_SVD::Solve");
        // OCCT L44: Standard_DimensionError_Raise_if((RowA != B.Length()) ||
        //            (X.Length() != Diag.Length()), " ");
        assert!(
            self.row_a == b.len() && x.len() == self.diag.len(),
            "Standard_DimensionError: math_SVD::Solve"
        );

        // OCCT L46-48: math_Vector BB(1, U.RowNumber()); BB.Init(0.0);
        //              BB.Set(1, B.Length(), B);
        let mut bb = VecD::new(self.u.n_rows());
        for i in 1..=bb.len() {
            bb.set(i, 0.0);
        }
        for i in 1..=b.len() {
            bb.set(i, b.get(i));
        }
        // OCCT L49: double wmin = Eps * Diag(Diag.Max()); — Max() is the
        // 1-based index of the maximum singular value.
        let wmin = eps * self.diag.get(self.diag_max_index());
        // OCCT L50-56.
        for i in 1..=self.diag.len() {
            if self.diag.get(i) < wmin {
                self.diag.set(i, 0.0);
            }
        }

        // OCCT L58-67: handle custom bounds in the X vector — VecD is always
        // 1-based (architecture difference), so X.Lower() == 1 always holds
        // and the direct SVD_Solve(U, Diag, V, BB, X) path applies.
        svd_solve(&self.u, &self.diag, &self.v, &bb, x);
    }

    /// OCCT math_SVD::PseudoInverse(Inv, Eps) (math_SVD.cxx L70-102) —
    /// computes the inverse Inv of matrix A such as A * Inv = Identity.
    pub fn pseudo_inverse(&mut self, result: &mut MatD, eps: f64) {
        // OCCT L74: StdFail_NotDone_Raise_if(!Done, " ");
        assert!(self.done, "StdFail_NotDone: math_SVD::PseudoInverse");

        // OCCT L76-83.
        let wmin = eps * self.diag.get(self.diag_max_index());
        for i in 1..=self.diag.len() {
            if self.diag.get(i) < wmin {
                self.diag.set(i, 0.0);
            }
        }

        // OCCT L85-87: int ColA = Diag.Length(); math_Vector VNorme(1,
        //              U.RowNumber()); math_Vector Column(1, ColA);
        let col_a = self.diag.len();
        let mut v_norme = VecD::new(self.u.n_rows());
        let mut column = VecD::new(col_a);

        // OCCT L89-101.
        for j in 1..=self.row_a {
            for i in 1..=v_norme.len() {
                v_norme.set(i, 0.0);
            }
            v_norme.set(j, 1.0);
            svd_solve(&self.u, &self.diag, &self.v, &v_norme, &mut column);
            for i in 1..=col_a {
                result.set(i, j, column.get(i));
            }
        }
    }

    /// OCCT math_SVD::Dump(o) (math_SVD.cxx L104-116) — prints information on
    /// the current state of the object. Standard_OStream has no Rust
    /// equivalent here, so the dump text is returned as a String
    /// (architecture difference).
    pub fn dump(&self) -> String {
        let mut text = String::from("math_SVD");
        if self.done {
            text.push_str(" Status = Done \n");
        } else {
            text.push_str(" Status = not Done \n");
        }
        text
    }
}
