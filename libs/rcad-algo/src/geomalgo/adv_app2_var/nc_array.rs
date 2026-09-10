//! Architecture carrier for `NCollection_Array2<T>` (OCCT TKFoundation
//! NCollection) inside the AdvApp2Var class layer.
//!
//! OCCT stores an Array2 with explicit lower/upper bounds per dimension and
//! addresses it as `Value(i, j)` with the caller's 1-based-or-lower indices.
//! rcad keeps the OCCT addressing (the logical indices are stored bounds, the
//! data is a flat row-major buffer); this is the encoding of the OCCT
//! `NCollection_Array2` / `NCollection_HArray2` used by AdvApp2Var_Node,
//! AdvApp2Var_Iso, AdvApp2Var_Patch and AdvApp2Var_ApproxAFunc2Var.
//!
//! `occ::handle<NCollection_HArray2<double>>` (a nullable shared handle) is
//! modeled as `Option<Array2<f64>>` at the member sites.

/// OCCT NCollection_Array2 - fixed-size 2D array with explicit bounds.
#[derive(Debug, Clone, PartialEq)]
pub struct Array2<T: Copy> {
    /// OCCT myLowerRow.
    row_lower: i32,
    /// OCCT myUpperRow.
    row_upper: i32,
    /// OCCT myLowerCol.
    col_lower: i32,
    /// OCCT myUpperCol.
    col_upper: i32,
    /// Flat row-major storage.
    data: Vec<T>,
}

impl<T: Copy + Default> Array2<T> {
    /// OCCT NCollection_Array2(theLowerRow, theUpperRow, theLowerCol,
    /// theUpperCol) - allocates without initializing (the OCCT ctor leaves
    /// the memory uninitialized; callers either Init() or SetValue() every
    /// entry, which rcad reproduces with Default placeholders).
    pub fn new(row_lower: i32, row_upper: i32, col_lower: i32, col_upper: i32) -> Self {
        let nrow = (row_upper - row_lower + 1).max(0) as usize;
        let ncol = (col_upper - col_lower + 1).max(0) as usize;
        Array2 {
            row_lower,
            row_upper,
            col_lower,
            col_upper,
            data: vec![T::default(); nrow * ncol],
        }
    }

    /// OCCT Init(theValue) - fills the whole array.
    pub fn init(&mut self, value: T) {
        for x in &mut self.data {
            *x = value;
        }
    }

    fn flat(&self, row: i32, col: i32) -> usize {
        // OCCT: (theRow - myLowerRow) * RowLength() + (theCol - myLowerCol).
        ((row - self.row_lower) * (self.col_upper - self.col_lower + 1)
            + (col - self.col_lower)) as usize
    }

    /// OCCT Value(theRow, theCol).
    pub fn value(&self, row: i32, col: i32) -> T {
        self.data[self.flat(row, col)]
    }

    /// OCCT SetValue(theRow, theCol, theValue).
    pub fn set_value(&mut self, row: i32, col: i32, value: T) {
        let k = self.flat(row, col);
        self.data[k] = value;
    }

    /// The flat row-major buffer (the `->Array2()` pass-through the OCCT
    /// handle wrappers expose to the Convert_GridPolynomialToPoles ctor).
    pub fn data(&self) -> &[T] {
        &self.data
    }

    /// OCCT ChangeValue(theRow, theCol) - mutable access.
    pub fn change_value(&mut self, row: i32, col: i32) -> &mut T {
        let k = self.flat(row, col);
        &mut self.data[k]
    }

    /// OCCT RowLength() - number of columns.
    pub fn row_length(&self) -> i32 {
        self.col_upper - self.col_lower + 1
    }

    /// OCCT ColLength() - number of rows.
    pub fn col_length(&self) -> i32 {
        self.row_upper - self.row_lower + 1
    }
}
