// OCCT Math* packages
pub mod adv_approx;
pub mod bnd;
pub mod bspl;
pub mod bspl_lib;
pub mod bvh;
pub mod convert_comp_polynomial_to_poles;
pub mod direct_polynomial_roots;
pub mod el;
pub mod gauss_points;
pub mod gauss_tables;
pub mod gp;
pub mod math_gauss;
pub mod math_poly;
pub mod newton_function_root;
pub mod p_lib_jacobi;
pub mod p_lib_jacobi_data;
pub mod plib;
pub mod root;
pub mod poly;
pub mod opt;
pub mod lin;
pub mod integ;
pub mod sys;

// keep flat — OCCT LProp package (partial)
pub mod top_loc;
pub mod curvature;

// OCCT CSLib package (surface normal computation)
pub mod cs_lib;

// legacy flat modules (keep during migration)
pub mod arc_length;
pub mod fit;
pub mod math_utils;
pub mod projection;
pub mod properties;

// =============================================================================
// math_Vector / math_Matrix / math_IntegerVector — 1:1 of the OCCT TKMath
// containers (1-based indexing).
// =============================================================================

/// OCCT math_Vector: a 1-based double array.
#[derive(Debug, Clone)]
pub struct VecD {
    pub v: Vec<f64>,
}

impl VecD {
    pub fn new(len: usize) -> Self {
        VecD { v: vec![0.0; len] }
    }
    /// 1-based index.
    pub fn get(&self, i: usize) -> f64 {
        self.v[i - 1]
    }
    pub fn set(&mut self, i: usize, x: f64) {
        self.v[i - 1] = x;
    }
    pub fn len(&self) -> usize {
        self.v.len()
    }
    /// math_Vector::Multiplied (dot product).
    pub fn multiplied(&self, other: &VecD) -> f64 {
        self.v.iter().zip(other.v.iter()).map(|(a, b)| a * b).sum()
    }
    /// math_Vector::Norm2.
    pub fn norm2(&self) -> f64 {
        self.v.iter().map(|x| x * x).sum()
    }
}

/// OCCT math_Matrix: a 1-based double matrix (rows then columns).
#[derive(Debug, Clone)]
pub struct MatD {
    pub m: Vec<Vec<f64>>,
}

impl MatD {
    pub fn new(n_rows: usize, n_cols: usize) -> Self {
        MatD {
            m: vec![vec![0.0; n_cols]; n_rows],
        }
    }
    pub fn get(&self, r: usize, c: usize) -> f64 {
        self.m[r - 1][c - 1]
    }
    pub fn set(&mut self, r: usize, c: usize, x: f64) {
        self.m[r - 1][c - 1] = x;
    }
    pub fn n_rows(&self) -> usize {
        self.m.len()
    }
    pub fn n_cols(&self) -> usize {
        self.m[0].len()
    }
}

/// OCCT math_IntegerVector: a 1-based int array.
#[derive(Debug, Clone)]
pub struct IntVec {
    pub v: Vec<i32>,
}

impl IntVec {
    pub fn new(len: usize) -> Self {
        IntVec { v: vec![0; len] }
    }
    pub fn get(&self, i: usize) -> i32 {
        self.v[i - 1]
    }
    pub fn set(&mut self, i: usize, x: i32) {
        self.v[i - 1] = x;
    }
    pub fn len(&self) -> usize {
        self.v.len()
    }
}

/// OCCT GeomAbs_Shape (TKG3d GeomAbs package) — continuity order.  Placed
/// here because the TKMath approximation stack (PLib / AdvApprox) consumes
/// it; re-exported for the algorithm crates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsShape {
    C0,
    C1,
    C2,
    C3,
    CN,
}

#[cfg(test)]
pub mod math_gtests;
#[cfg(test)]
pub mod tkmath_gtests;
