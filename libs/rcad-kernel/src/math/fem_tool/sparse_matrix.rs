//! OCCT FEmTool_SparseMatrix (ModelingData/TKGeomBase/FEmTool).
//!
//! 1:1 translation of `FEmTool_SparseMatrix.hxx` (L17-66): the abstract
//! sparse-matrix interface of the FEmTool package.  The OCCT
//! `Standard_Transient` base becomes a Rust trait (abstract bases map to
//! traits per the kernel convention); the single implementation in the
//! package is [`super::ProfileMatrix`].

use super::super::math_matrix::Vector;

/// OCCT FEmTool_SparseMatrix - sparse matrix definition
/// (FEmTool_SparseMatrix.hxx L28-64).
pub trait SparseMatrix {
    /// OCCT FEmTool_SparseMatrix::Init(Value).
    fn init(&mut self, value: f64);

    /// OCCT FEmTool_SparseMatrix::ChangeValue(I, J) - returns `double&`.
    fn change_value(&mut self, i: i32, j: i32) -> &mut f64;

    /// OCCT FEmTool_SparseMatrix::Decompose - to make a factorization of
    /// `self`.
    fn decompose(&mut self) -> bool;

    /// OCCT FEmTool_SparseMatrix::Solve(B, X) - direct solve of A X = B.
    fn solve(&self, b: &Vector, x: &mut Vector);

    /// OCCT FEmTool_SparseMatrix::Prepare - make preparation to iterative
    /// solve.
    fn prepare(&mut self) -> bool;

    /// OCCT FEmTool_SparseMatrix::Solve(B, Init, X, Residual, Tolerance,
    /// NbIterations) - iterative solve of A X = B.
    fn solve_iterative(
        &self,
        b: &Vector,
        init: &Vector,
        x: &mut Vector,
        residual: &mut Vector,
        tolerance: f64,
        nb_iterations: i32,
    );

    /// OCCT FEmTool_SparseMatrix::Multiplied(X, MX) - returns the product of
    /// a SparseMatrix by a vector.  An exception is raised if the dimensions
    /// are different.
    fn multiplied(&self, x: &Vector, mx: &mut Vector);

    /// OCCT FEmTool_SparseMatrix::RowNumber - returns the row range of a
    /// matrix.
    fn row_number(&self) -> i32;

    /// OCCT FEmTool_SparseMatrix::ColNumber - returns the column range of
    /// the matrix.
    fn col_number(&self) -> i32;
}
