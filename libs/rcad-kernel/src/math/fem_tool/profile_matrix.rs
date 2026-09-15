//! OCCT FEmTool_ProfileMatrix (ModelingData/TKGeomBase/FEmTool).
//!
//! 1:1 translation of `FEmTool_ProfileMatrix.hxx` (L17-85) and
//! `FEmTool_ProfileMatrix.cxx` (L1-322).
//!
//! Container mapping: the OCCT flat 1-based storages (`ProfileMatrix`,
//! `SMatrix`, `NextCoeff` HArray1 with the `SMA--` / `PM--` / `NC--`
//! pointer trick) become `Vec`s indexed with explicit `k - 1`; the
//! `profile` NCollection_Array2<int>(1, 2, 1, n) becomes a
//! `Vec<[i32; 2]>` (`profile(1, i)` is `.0`, `profile(2, i)` is `.1`).
//! Vectors keep the bound-aware `math_matrix::Vector` (math_Vector).

use super::super::math_matrix::Vector;
use super::sparse_matrix::SparseMatrix;

/// OCCT FEmTool_ProfileMatrix - symmetric sparse profile matrix useful for
/// 1D finite element methods (FEmTool_ProfileMatrix.hxx L32-83).
#[derive(Debug, Clone)]
pub struct ProfileMatrix {
    /// OCCT `profile`: NCollection_Array2<int>(1, 2, 1, FirstIndexes.Length()).
    /// `profile(1, i)` = half-bandwidth of row i; `profile(2, i)` = end
    /// address of row i in the flat storages.
    profile: Vec<[i32; 2]>,
    /// OCCT `ProfileMatrix` HArray1<double>(1, profile(2, n)).
    profile_matrix: Vec<f64>,
    /// OCCT `SMatrix` HArray1<double>(1, profile(2, n)).
    s_matrix: Vec<f64>,
    /// OCCT `NextCoeff` HArray1<int>(1, profile(2, n)).
    next_coeff: Vec<i32>,
    /// OCCT `IsDecomp`.
    is_decomp: bool,
}

impl ProfileMatrix {
    /// OCCT FEmTool_ProfileMatrix::FEmTool_ProfileMatrix
    /// (FEmTool_ProfileMatrix.cxx L32-69).  `first_indexes` is the OCCT
    /// NCollection_Array1<int>(1, n) `FirstIndexes`; the slice is 0-based,
    /// so OCCT `FirstIndexes(i)` is `first_indexes[i - 1]`.
    pub fn new(first_indexes: &[i32]) -> Self {
        let n = first_indexes.len() as i32;

        let mut profile = vec![[0i32; 2]; n as usize];
        // profile(1, 1) = 0; profile(2, 1) = 1;
        profile[0][0] = 0;
        profile[0][1] = 1;
        for i in 2..=n {
            // profile(1, i) = i - FirstIndexes(i);
            profile[(i - 1) as usize][0] = i - first_indexes[(i - 1) as usize];
            // profile(2, i) = profile(2, i - 1) + profile(1, i) + 1;
            profile[(i - 1) as usize][1] = profile[(i - 2) as usize][1] + profile[(i - 1) as usize][0] + 1;
        }
        let total = profile[(n - 1) as usize][1] as usize;

        let mut next_coeff = vec![0i32; total];

        let mut k = 0usize;
        for i in 1..=n {
            let mut j = first_indexes[(i - 1) as usize];
            while j <= i {
                let mut l = i + 1;
                while l <= n && j < first_indexes[(l - 1) as usize] {
                    l += 1;
                }

                if l > n {
                    next_coeff[k] = 0;
                } else {
                    next_coeff[k] = l;
                }
                k += 1;
                j += 1;
            }
        }

        ProfileMatrix {
            profile,
            profile_matrix: vec![0.0; total],
            s_matrix: vec![0.0; total],
            next_coeff,
            is_decomp: false,
        }
    }

    /// OCCT FEmTool_ProfileMatrix::IsInProfile (cxx L264-276).
    pub fn is_in_profile(&self, i: i32, j: i32) -> bool {
        if j <= i {
            (i - j) <= self.profile[(i - 1) as usize][0]
        } else if (j - i) <= self.profile[(j - 1) as usize][0] {
            true
        } else {
            false
        }
    }

    /// OCCT FEmTool_ProfileMatrix::OutM (cxx L278-302) - debug print of the
    /// matrix A and the NextCoeff chain.
    pub fn out_m(&self) {
        println!("Matrix A");
        for i in 1..=self.row_number() {
            for _j in 1..(i - self.profile[(i - 1) as usize][0]) {
                print!("0 ");
            }

            for j in (self.profile[(i - 1) as usize][1] - self.profile[(i - 1) as usize][0])
                ..=self.profile[(i - 1) as usize][1]
            {
                print!("{} ", self.profile_matrix[(j - 1) as usize]);
            }
            println!();
        }

        println!("NextCoeff");
        for i in 1..=self.profile[(self.row_number() - 1) as usize][1] {
            print!("{} ", self.next_coeff[(i - 1) as usize]);
        }
        println!();
    }

    /// OCCT FEmTool_ProfileMatrix::OutS (cxx L304-321) - debug print of the
    /// decomposed matrix S.
    pub fn out_s(&self) {
        println!("Matrix S");
        for i in 1..=self.row_number() {
            for _j in 1..(i - self.profile[(i - 1) as usize][0]) {
                print!("0 ");
            }

            for j in (self.profile[(i - 1) as usize][1] - self.profile[(i - 1) as usize][0])
                ..=self.profile[(i - 1) as usize][1]
            {
                print!("{} ", self.s_matrix[(j - 1) as usize]);
            }
            println!();
        }
    }
}

impl SparseMatrix for ProfileMatrix {
    /// OCCT FEmTool_ProfileMatrix::Init (cxx L73-77).
    fn init(&mut self, value: f64) {
        self.profile_matrix.fill(value);
        self.is_decomp = false;
    }

    /// OCCT FEmTool_ProfileMatrix::ChangeValue (cxx L81-97).
    fn change_value(&mut self, i: i32, j: i32) -> &mut f64 {
        let mut ind = i - j;
        if ind < 0 {
            ind = -ind;
            assert!(
                ind <= self.profile[(j - 1) as usize][0],
                "Standard_OutOfRange: FEmTool_ProfileMatrix::ChangeValue"
            );
            ind = self.profile[(j - 1) as usize][1] - ind;
        } else {
            assert!(
                ind <= self.profile[(i - 1) as usize][0],
                "Standard_OutOfRange: FEmTool_ProfileMatrix::ChangeValue"
            );
            ind = self.profile[(i - 1) as usize][1] - ind;
        }
        &mut self.profile_matrix[(ind - 1) as usize]
    }

    /// OCCT FEmTool_ProfileMatrix::Decompose (cxx L101-150) - Cholesky-like
    /// factorization A = S t(S); returns false if the matrix is not positive
    /// defined.
    fn decompose(&mut self) -> bool {
        let eps = 1.0e-32;

        // SMatrix->Init(0.); the OCCT `SMA--` / `PM--` pointer offsets are
        // translated as explicit (k - 1) indexing.
        self.s_matrix.fill(0.0);
        self.is_decomp = false;
        for j in 1..=self.row_number() {
            let diag_addr = self.profile[(j - 1) as usize][1];
            let kj = j - self.profile[(j - 1) as usize][0];
            let mut sum = 0.0;
            let mut k = diag_addr - self.profile[(j - 1) as usize][0];
            while k < diag_addr {
                // Sum += SMA[k] * SMA[k];
                sum += self.s_matrix[(k - 1) as usize] * self.s_matrix[(k - 1) as usize];
                k += 1;
            }

            let mut a = self.profile_matrix[(diag_addr - 1) as usize] - sum;
            if a < eps {
                return false; // Matrix is not positive defined
            }
            a = a.sqrt();
            self.s_matrix[(diag_addr - 1) as usize] = a;

            let mut curr_addr = diag_addr;
            loop {
                // while ((i = NextCoeff->Value(CurrAddr)) > 0)
                let i = self.next_coeff[(curr_addr - 1) as usize];
                if i <= 0 {
                    break;
                }
                curr_addr = self.profile[(i - 1) as usize][1] - (i - j);

                // Computation of Sum of S_ik . S_jk for k = 1,..,j-1
                let mut sum = 0.0;
                let kmin = (i - self.profile[(i - 1) as usize][0]).max(kj);
                let mut ik = self.profile[(i - 1) as usize][1] - i + kmin;
                let mut jk = diag_addr - j + kmin;
                let mut k = kmin;
                while k < j {
                    sum += self.s_matrix[(ik - 1) as usize] * self.s_matrix[(jk - 1) as usize];
                    ik += 1;
                    jk += 1;
                    k += 1;
                }
                self.s_matrix[(curr_addr - 1) as usize] =
                    (self.profile_matrix[(curr_addr - 1) as usize] - sum) / a;
            }
        }
        self.is_decomp = true;
        self.is_decomp
    }

    /// OCCT FEmTool_ProfileMatrix::Solve (cxx L156-201) - resolution of the
    /// system S t(S) X = B.
    fn solve(&self, b: &Vector, x: &mut Vector) {
        assert!(
            self.is_decomp,
            "StdFail_NotDone: Decomposition must be done"
        );

        // Resolution of Sw = B;
        for i in 1..=self.row_number() {
            let diag_addr = self.profile[(i - 1) as usize][1];
            let mut sum = 0.0;
            let mut j = i - self.profile[(i - 1) as usize][0];
            let mut jj = diag_addr - (i - j);
            while j < i {
                sum += self.s_matrix[(jj - 1) as usize] * x.get(j);
                j += 1;
                jj += 1;
            }
            x.set(i, (b.get(i) - sum) / self.s_matrix[(diag_addr - 1) as usize]);
        }

        // Resolution of t(S)X = w;
        for i in (1..=self.col_number()).rev() {
            let diag_addr = self.profile[(i - 1) as usize][1];
            let mut j = self.next_coeff[(diag_addr - 1) as usize];
            let mut sum = 0.0;
            while j > 0 {
                let curr_addr = self.profile[(j - 1) as usize][1] - (j - i);
                sum += self.s_matrix[(curr_addr - 1) as usize] * x.get(j);
                j = self.next_coeff[(curr_addr - 1) as usize];
            }
            x.set(i, (x.get(i) - sum) / self.s_matrix[(diag_addr - 1) as usize]);
        }
    }

    /// OCCT FEmTool_ProfileMatrix::Prepare (cxx L203-206).
    fn prepare(&mut self) -> bool {
        panic!("Standard_NotImplemented: FEmTool_ProfileMatrix::Prepare");
    }

    /// OCCT FEmTool_ProfileMatrix::Solve (iterative, cxx L210-218).
    fn solve_iterative(
        &self,
        _b: &Vector,
        _init: &Vector,
        _x: &mut Vector,
        _residual: &mut Vector,
        _tolerance: f64,
        _nb_iterations: i32,
    ) {
        panic!("Standard_NotImplemented: FEmTool_ProfileMatrix::Solve");
    }

    /// OCCT FEmTool_ProfileMatrix::Multiplied (cxx L224-252) - MX = H*X.
    fn multiplied(&self, x: &Vector, mx: &mut Vector) {
        for i in 1..=self.row_number() {
            let diag_addr = self.profile[(i - 1) as usize][1];
            let mut m_i = 0.0;
            let mut j = i - self.profile[(i - 1) as usize][0];
            let mut jj = diag_addr - (i - j);
            while j <= i {
                m_i += self.profile_matrix[(jj - 1) as usize] * x.get(j);
                j += 1;
                jj += 1;
            }

            let mut curr_addr = diag_addr;
            let mut j = self.next_coeff[(curr_addr - 1) as usize];
            while j > 0 {
                curr_addr = self.profile[(j - 1) as usize][1] - (j - i);
                m_i += self.profile_matrix[(curr_addr - 1) as usize] * x.get(j);
                j = self.next_coeff[(curr_addr - 1) as usize];
            }
            mx.set(i, m_i);
        }
    }

    /// OCCT FEmTool_ProfileMatrix::RowNumber (cxx L254-257).
    fn row_number(&self) -> i32 {
        self.profile.len() as i32
    }

    /// OCCT FEmTool_ProfileMatrix::ColNumber (cxx L259-262).
    fn col_number(&self) -> i32 {
        self.profile.len() as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ProfileMatrix solve, hand-solved for the 3x3 system
    ///   A = [[2,1,0],[1,3,1],[0,1,2]],  b = [3,5,3]  =>  x = [1,1,1].
    /// Cholesky hand-check: S = [[sqrt(2),0,0],[1/sqrt(2),sqrt(5/2),0],
    /// [0,1/sqrt(5/2),sqrt(8/5)]]; forward/backward substitution gives
    /// exactly [1,1,1].  FirstIndexes = [1,1,2] (row 3 starts at column 2).
    #[test]
    fn profile_matrix_3x3_hand_solved() {
        // OCCT NCollection_Array1<int> FirstIndexes(1, 3) = {1, 1, 2}.
        let mut a = ProfileMatrix::new(&[1, 1, 2]);
        *a.change_value(1, 1) = 2.0;
        *a.change_value(2, 1) = 1.0;
        *a.change_value(2, 2) = 3.0;
        *a.change_value(3, 2) = 1.0;
        *a.change_value(3, 3) = 2.0;

        // Symmetric addressing: ChangeValue(2, 3) is the same cell as
        // ChangeValue(3, 2).
        assert!((*a.change_value(2, 3) - 1.0).abs() < f64::EPSILON);

        assert!(a.decompose(), "matrix is positive defined");

        let mut b = Vector::new_init(1, 3, 0.0);
        b.set(1, 3.0);
        b.set(2, 5.0);
        b.set(3, 3.0);
        let mut x = Vector::new(1, 3);
        a.solve(&b, &mut x);

        assert!((x.get(1) - 1.0).abs() < 1e-13, "x1 = {}", x.get(1));
        assert!((x.get(2) - 1.0).abs() < 1e-13, "x2 = {}", x.get(2));
        assert!((x.get(3) - 1.0).abs() < 1e-13, "x3 = {}", x.get(3));

        // Multiplied: A*x must reproduce b.
        let mut mx = Vector::new(1, 3);
        a.multiplied(&x, &mut mx);
        assert!((mx.get(1) - 3.0).abs() < 1e-13);
        assert!((mx.get(2) - 5.0).abs() < 1e-13);
        assert!((mx.get(3) - 3.0).abs() < 1e-13);
    }
}
