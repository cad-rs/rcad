//! OCCT MAT2d_BiInt — a set of two integers.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_BiInt.hxx L17-66, MAT2d_BiInt.cxx L19-48
//!
//! Used as the key of MAT2d_Circuit::linkRefEqui (NCollection_DataMap
//! keyed on MAT2d_BiInt), hence Hash + Eq mirroring the std::hash
//! specialization of the .hxx.

/// OCCT MAT2d_BiInt (MAT2d_BiInt.hxx L28-50) — a set of two integers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Mat2dBiInt {
    i1: i32,
    i2: i32,
}

impl Mat2dBiInt {
    /// OCCT MAT2d_BiInt::MAT2d_BiInt(I1, I2) (cxx L19-23).
    pub fn new(i1: i32, i2: i32) -> Self {
        Mat2dBiInt { i1, i2 }
    }

    /// OCCT MAT2d_BiInt::FirstIndex() const (cxx L25-28).
    pub fn first_index(&self) -> i32 {
        self.i1
    }

    /// OCCT MAT2d_BiInt::SecondIndex() const (cxx L30-33).
    pub fn second_index(&self) -> i32 {
        self.i2
    }

    /// OCCT MAT2d_BiInt::FirstIndex(I1) (cxx L35-38).
    pub fn set_first_index(&mut self, i1: i32) {
        self.i1 = i1;
    }

    /// OCCT MAT2d_BiInt::SecondIndex(I2) (cxx L40-43).
    pub fn set_second_index(&mut self, i2: i32) {
        self.i2 = i2;
    }

    /// OCCT MAT2d_BiInt::IsEqual(B) const (cxx L45-48).
    /// (operator== of the .hxx L45 delegates here; PartialEq mirrors it.)
    pub fn is_equal(&self, b: &Mat2dBiInt) -> bool {
        self.i1 == b.first_index() && self.i2 == b.second_index()
    }
}
