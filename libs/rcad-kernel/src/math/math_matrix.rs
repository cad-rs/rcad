// OCCT math_Matrix / math_Vector / math_IntegerVector (TKMath) — 1-based
// storage with arbitrary LowerRow/LowerCol bounds. The shells MatD/VecD/
// IntVec in math/mod.rs are the bounds-less legacy forms; these are the
// bound-aware forms required by math_Crout / math_Uzawa and the
// AppParCurves approximation templates.

use super::{IntVec, MatD, VecD};

/// OCCT math_Matrix.
#[derive(Debug, Clone)]
pub struct Matrix {
    pub data: MatD,
    pub lower_row: i32,
    pub lower_col: i32,
}

impl Matrix {
    /// OCCT math_Matrix(I1, I2, J1, J2) (zero-initialized storage).
    pub fn new(i1: i32, i2: i32, j1: i32, j2: i32) -> Self {
        assert!(i2 >= i1 && j2 >= j1, "math_Matrix: bad range");
        Matrix {
            data: MatD::new((i2 - i1 + 1) as usize, (j2 - j1 + 1) as usize),
            lower_row: i1,
            lower_col: j1,
        }
    }

    /// OCCT math_Matrix(I1, I2, J1, J2, InitValue).
    pub fn new_init(i1: i32, i2: i32, j1: i32, j2: i32, init: f64) -> Self {
        let mut m = Matrix::new(i1, i2, j1, j2);
        m.init(init);
        m
    }

    /// OCCT Init(Value).
    pub fn init(&mut self, init: f64) {
        for r in 1..=self.row_number() {
            for c in 1..=self.col_number() {
                self.data.m[(r - 1) as usize][(c - 1) as usize] = init;
            }
        }
    }

    /// OCCT math_Matrix::Set(I1, I2, J1, J2, M) — replaces the sub-matrix
    /// [I1..I2, J1..J2] by the contents of M.
    pub fn set_block(&mut self, i1: i32, i2: i32, j1: i32, j2: i32, m: &Matrix) {
        for i in i1..=i2 {
            for j in j1..=j2 {
                let v = m.get(i - i1 + m.lower_row, j - j1 + m.lower_col);
                self.set(i, j, v);
            }
        }
    }

    #[inline]
    pub fn get(&self, i: i32, j: i32) -> f64 {
        self.data.get((i - self.lower_row + 1) as usize, (j - self.lower_col + 1) as usize)
    }

    #[inline]
    pub fn set(&mut self, i: i32, j: i32, v: f64) {
        self.data
            .set((i - self.lower_row + 1) as usize, (j - self.lower_col + 1) as usize, v);
    }

    /// OCCT RowNumber().
    pub fn row_number(&self) -> i32 {
        self.data.n_rows() as i32
    }

    /// OCCT ColNumber().
    pub fn col_number(&self) -> i32 {
        self.data.n_cols() as i32
    }

    /// OCCT LowerRow().
    pub fn lower_row(&self) -> i32 {
        self.lower_row
    }

    /// OCCT UpperRow().
    pub fn upper_row(&self) -> i32 {
        self.lower_row + self.row_number() - 1
    }

    /// OCCT LowerCol().
    pub fn lower_col(&self) -> i32 {
        self.lower_col
    }

    /// OCCT UpperCol().
    pub fn upper_col(&self) -> i32 {
        self.lower_col + self.col_number() - 1
    }

    /// Normalized (1-based) storage view, for the legacy MatD-based entry
    /// points (e.g. rcad Householder).
    pub fn data(&self) -> &MatD {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut MatD {
        &mut self.data
    }
}

/// OCCT math_Vector.
#[derive(Debug, Clone)]
pub struct Vector {
    pub data: VecD,
    pub lower: i32,
}

impl Vector {
    /// OCCT math_Vector(I1, I2).
    pub fn new(i1: i32, i2: i32) -> Self {
        assert!(i2 >= i1, "math_Vector: bad range");
        Vector {
            data: VecD::new((i2 - i1 + 1) as usize),
            lower: i1,
        }
    }

    /// OCCT math_Vector(I1, I2, InitValue).
    pub fn new_init(i1: i32, i2: i32, init: f64) -> Self {
        let mut v = Vector::new(i1, i2);
        for r in 1..=v.data.len() {
            v.data.set(r, init);
        }
        v
    }

    #[inline]
    pub fn get(&self, i: i32) -> f64 {
        self.data.get((i - self.lower + 1) as usize)
    }

    #[inline]
    pub fn set(&mut self, i: i32, v: f64) {
        self.data.set((i - self.lower + 1) as usize, v);
    }

    /// OCCT Length().
    pub fn length(&self) -> i32 {
        self.data.len() as i32
    }

    /// OCCT Lower().
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// OCCT Upper().
    pub fn upper(&self) -> i32 {
        self.lower + self.length() - 1
    }
}

/// OCCT math_IntegerVector.
#[derive(Debug, Clone)]
pub struct IntegerVector {
    pub data: IntVec,
    pub lower: i32,
}

impl IntegerVector {
    /// OCCT math_IntegerVector(I1, I2, InitValue).
    pub fn new_init(i1: i32, i2: i32, init: i32) -> Self {
        assert!(i2 >= i1, "math_IntegerVector: bad range");
        let mut v = IntegerVector {
            data: IntVec::new((i2 - i1 + 1) as usize),
            lower: i1,
        };
        for r in 1..=v.data.len() {
            v.data.set(r, init);
        }
        v
    }

    #[inline]
    pub fn get(&self, i: i32) -> i32 {
        self.data.get((i - self.lower + 1) as usize)
    }

    #[inline]
    pub fn set(&mut self, i: i32, v: i32) {
        self.data.set((i - self.lower + 1) as usize, v);
    }

    /// OCCT Length().
    pub fn length(&self) -> i32 {
        self.data.len() as i32
    }

    /// OCCT Lower().
    pub fn lower(&self) -> i32 {
        self.lower
    }
}
