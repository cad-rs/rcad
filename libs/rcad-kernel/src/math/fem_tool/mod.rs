//! OCCT FEmTool package (ModelingData/TKGeomBase/FEmTool) - finite element
//! tooling for variational approximation: sparse profile solver, Hermit-
//! Jacobi element curves, system assembly and elementary (linear) criteria.
//!
//! Package classes, each in its own child module with OCCT line anchors:
//! `FEmTool_SparseMatrix` (sparse_matrix), `FEmTool_ProfileMatrix`
//! (profile_matrix), `FEmTool_Curve` (curve), `FEmTool_Assembly`
//! (assembly), `FEmTool_ElementaryCriterion` (elementary_criterion),
//! `FEmTool_LinearTension` (linear_tension), `FEmTool_LinearFlexion`
//! (linear_flexion), `FEmTool_LinearJerk` (linear_jerk),
//! `FEmTool_ElementsOfRefMatrix` (elements_of_ref_matrix).
//!
//! Hosted companions (documented in the child modules): `math_FunctionSet`
//! + `math_GaussSetIntegration` (gauss_set_integration), required by the
//! criteria constructors; the `math` module registration is closed so they
//! travel with their only consumers.  Also hosted here are the minimal
//! NCollection mirrors used by the package: `IntArray1` /
//! `IntArray2` (int arrays with arbitrary bounds) and `AssemblingTable`
//! (`HArray2<Handle(HArray1<int>)>`).

pub mod assembly;
pub mod curve;
pub mod elementary_criterion;
pub mod elements_of_ref_matrix;
pub mod gauss_set_integration;
pub mod linear_flexion;
pub mod linear_jerk;
pub mod linear_tension;
pub mod profile_matrix;
pub mod sparse_matrix;

pub use assembly::Assembly;
pub use curve::Curve;
pub use elementary_criterion::{CoeffHandle, CriterionData, ElementaryCriterion};
pub use elements_of_ref_matrix::ElementsOfRefMatrix;
pub use gauss_set_integration::{FunctionSet, GaussSetIntegration};
pub use linear_flexion::LinearFlexion;
pub use linear_jerk::LinearJerk;
pub use linear_tension::LinearTension;
pub use profile_matrix::ProfileMatrix;
pub use sparse_matrix::SparseMatrix;

/// Mirror of `NCollection_Array1<int>` / `NCollection_HArray1<int>`: a
/// 1-based-storage int array with an arbitrary lower bound (the kernel has
/// no generic NCollection port; only the int forms used by FEmTool are
/// mirrored here).
#[derive(Debug, Clone)]
pub struct IntArray1 {
    lower: i32,
    data: Vec<i32>,
}

impl IntArray1 {
    /// OCCT NCollection_HArray1<int>(Lower, Upper).
    pub fn new(lower: i32, upper: i32) -> Self {
        assert!(upper >= lower - 1, "NCollection_HArray1: bad range");
        IntArray1 {
            lower,
            data: vec![0i32; (upper - lower + 1).max(0) as usize],
        }
    }

    /// OCCT Init(Value).
    pub fn init(&mut self, value: i32) {
        self.data.fill(value);
    }

    /// OCCT Lower().
    #[inline]
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// OCCT Upper().
    #[inline]
    pub fn upper(&self) -> i32 {
        self.lower + self.data.len() as i32 - 1
    }

    /// OCCT Value(Index).
    #[inline]
    pub fn value(&self, index: i32) -> i32 {
        self.data[(index - self.lower) as usize]
    }

    /// OCCT SetValue(Index, Value).
    #[inline]
    pub fn set_value(&mut self, index: i32, value: i32) {
        self.data[(index - self.lower) as usize] = value;
    }
}

/// Mirror of `NCollection_Array2<int>` / `NCollection_HArray2<int>`: an int
/// matrix with arbitrary row/column bounds.  OCCT semantics
/// (NCollection_Array2.hxx L198-201): RowLength() is the number of columns
/// and ColLength() the number of rows.
#[derive(Debug, Clone)]
pub struct IntArray2 {
    lower_row: i32,
    lower_col: i32,
    nb_rows: usize,
    nb_columns: usize,
    data: Vec<i32>,
}

impl IntArray2 {
    /// OCCT NCollection_HArray2<int>(LowerRow, UpperRow, LowerCol, UpperCol)
    /// (zero-initialized).
    pub fn new(lower_row: i32, upper_row: i32, lower_col: i32, upper_col: i32) -> Self {
        assert!(upper_row >= lower_row && upper_col >= lower_col, "bad range");
        let nb_rows = (upper_row - lower_row + 1) as usize;
        let nb_columns = (upper_col - lower_col + 1) as usize;
        IntArray2 {
            lower_row,
            lower_col,
            nb_rows,
            nb_columns,
            data: vec![0i32; nb_rows * nb_columns],
        }
    }

    /// OCCT NCollection_HArray2<int>(..., InitValue).
    pub fn new_init(lower_row: i32, upper_row: i32, lower_col: i32, upper_col: i32, init: i32) -> Self {
        let mut a = IntArray2::new(lower_row, upper_row, lower_col, upper_col);
        a.data.fill(init);
        a
    }

    /// OCCT Value(I, J).
    #[inline]
    pub fn value(&self, i: i32, j: i32) -> i32 {
        self.data[((i - self.lower_row) as usize) * self.nb_columns + (j - self.lower_col) as usize]
    }

    /// OCCT SetValue(I, J, Value).
    #[inline]
    pub fn set_value(&mut self, i: i32, j: i32, value: i32) {
        self.data[((i - self.lower_row) as usize) * self.nb_columns + (j - self.lower_col) as usize] =
            value;
    }

    /// OCCT LowerRow().
    #[inline]
    pub fn lower_row(&self) -> i32 {
        self.lower_row
    }

    /// OCCT UpperRow().
    #[inline]
    pub fn upper_row(&self) -> i32 {
        self.lower_row + self.nb_rows as i32 - 1
    }

    /// OCCT LowerCol().
    #[inline]
    pub fn lower_col(&self) -> i32 {
        self.lower_col
    }

    /// OCCT UpperCol().
    #[inline]
    pub fn upper_col(&self) -> i32 {
        self.lower_col + self.nb_columns as i32 - 1
    }

    /// OCCT ColLength() - number of rows (NCollection_Array2.hxx L201).
    #[inline]
    pub fn col_length(&self) -> i32 {
        self.nb_rows as i32
    }

    /// OCCT RowLength() - number of columns (NCollection_Array2.hxx L198).
    #[inline]
    pub fn row_length(&self) -> i32 {
        self.nb_columns as i32
    }
}

/// Mirror of `NCollection_HArray2<Handle(NCollection_HArray1<int>)>` - the
/// FEmTool assembly table mapping (dimension, element) to the element-local
/// index arrays.
#[derive(Debug, Clone)]
pub struct AssemblingTable {
    lower_row: i32,
    lower_col: i32,
    nb_rows: usize,
    nb_columns: usize,
    cells: Vec<IntArray1>,
}

impl AssemblingTable {
    /// OCCT construction from per-cell arrays with a shared shape
    /// (LowerRow, UpperRow, LowerCol, UpperCol); each cell array gets the
    /// same (lower, upper) bounds given here.
    pub fn new(
        lower_row: i32,
        upper_row: i32,
        lower_col: i32,
        upper_col: i32,
        cell_lower: i32,
        cell_upper: i32,
    ) -> Self {
        assert!(upper_row >= lower_row && upper_col >= lower_col, "bad range");
        let nb_rows = (upper_row - lower_row + 1) as usize;
        let nb_columns = (upper_col - lower_col + 1) as usize;
        AssemblingTable {
            lower_row,
            lower_col,
            nb_rows,
            nb_columns,
            cells: (0..nb_rows * nb_columns)
                .map(|_| IntArray1::new(cell_lower, cell_upper))
                .collect(),
        }
    }

    /// OCCT Table->Value(Dim, El) (shared cell access).
    #[inline]
    pub fn value(&self, dim: i32, el: i32) -> &IntArray1 {
        &self.cells[((dim - self.lower_row) as usize) * self.nb_columns + (el - self.lower_col) as usize]
    }

    /// OCCT Table->ChangeValue(Dim, El).
    #[inline]
    pub fn change_value(&mut self, dim: i32, el: i32) -> &mut IntArray1 {
        &mut self.cells
            [((dim - self.lower_row) as usize) * self.nb_columns + (el - self.lower_col) as usize]
    }

    /// OCCT Table->LowerRow().
    #[inline]
    pub fn lower_row(&self) -> i32 {
        self.lower_row
    }

    /// OCCT Table->UpperRow().
    #[inline]
    pub fn upper_row(&self) -> i32 {
        self.lower_row + self.nb_rows as i32 - 1
    }

    /// OCCT Table->LowerCol().
    #[inline]
    pub fn lower_col(&self) -> i32 {
        self.lower_col
    }

    /// OCCT Table->UpperCol().
    #[inline]
    pub fn upper_col(&self) -> i32 {
        self.lower_col + self.nb_columns as i32 - 1
    }
}
