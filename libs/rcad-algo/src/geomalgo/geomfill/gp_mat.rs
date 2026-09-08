//! OCCT gp_Mat (TKMath/gp) — minimal 1:1 port of the members consumed by the
//! GeomFill convertors (QuasiAngular / Polynomial): SetRotation, SetCross,
//! Powered, Multiplied (scalar / matrix / column-vector) and the gp_XYZ
//! row-vector Multiply.
//!
//! OCCT gp_Mat stores `myMat[3][3]` where `matrix[i][j]` is `Value(i+1, j+1)`;
//! this port keeps the identical row-major layout.

use glam::DVec3;

/// OCCT gp_Mat.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpMat {
    /// `mat[i][j]` == OCCT `myMat[i][j]` == `Value(i+1, j+1)` (0-based here).
    pub mat: [[f64; 3]; 3],
}

impl GpMat {
    /// OCCT gp_Mat() — identity (gp_Mat.hxx L129: SetIdentity in the default
    /// constructor... the empty constructor leaves values uninitialized; the
    /// convertors always overwrite them, so identity is a safe neutral form).
    pub fn identity() -> Self {
        GpMat {
            mat: [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        }
    }

    /// OCCT gp_Mat(a11,a12,a13,a21,a22,a23,a31,a32,a33) row-major constructor.
    #[allow(clippy::too_many_arguments)]
    pub fn from_rows(
        a11: f64,
        a12: f64,
        a13: f64,
        a21: f64,
        a22: f64,
        a23: f64,
        a31: f64,
        a32: f64,
        a33: f64,
    ) -> Self {
        GpMat {
            mat: [
                [a11, a12, a13],
                [a21, a22, a23],
                [a31, a32, a33],
            ],
        }
    }

    /// OCCT gp_Mat::SetRotation (gp_Mat.cxx): Rodrigues' rotation formula
    /// R = I + sin(theta) K + (1 - cos(theta)) K^2 with K the skew matrix of
    /// the normalized axis.
    pub fn set_rotation(&mut self, axis: DVec3, ang: f64) {
        // OCCT: const gp_XYZ aV = theAxis.Normalized();
        let a_v = axis.normalize_or_zero();
        let a = a_v.x;
        let b = a_v.y;
        let c = a_v.z;
        let a_cos = ang.cos();
        let a_sin = ang.sin();
        let a_om_cos = 1.0 - a_cos; // One minus cosine
        let a2 = a * a;
        let b2 = b * b;
        let c2 = c * c;
        let ab = a * b;
        let ac = a * c;
        let bc = b * c;
        self.mat[0][0] = 1.0 + a_om_cos * (-(b2 + c2));
        self.mat[0][1] = a_om_cos * ab - a_sin * c;
        self.mat[0][2] = a_om_cos * ac + a_sin * b;
        self.mat[1][0] = a_om_cos * ab + a_sin * c;
        self.mat[1][1] = 1.0 + a_om_cos * (-(a2 + c2));
        self.mat[1][2] = a_om_cos * bc - a_sin * a;
        self.mat[2][0] = a_om_cos * ac - a_sin * b;
        self.mat[2][1] = a_om_cos * bc + a_sin * a;
        self.mat[2][2] = 1.0 + a_om_cos * (-(a2 + b2));
    }

    /// OCCT gp_Mat::SetCross (gp_Mat.cxx) — the cross-product operator of
    /// theRef: this * v == theRef ^ v.
    pub fn set_cross(&mut self, the_ref: DVec3) {
        let x = the_ref.x;
        let y = the_ref.y;
        let z = the_ref.z;
        self.mat[0][0] = 0.0;
        self.mat[1][1] = 0.0;
        self.mat[2][2] = 0.0;
        self.mat[0][1] = -z;
        self.mat[0][2] = y;
        self.mat[1][2] = -x;
        self.mat[1][0] = z;
        self.mat[2][0] = -y;
        self.mat[2][1] = x;
    }

    /// OCCT gp_Mat::Powered(N) — matrix to the power N (gp_Mat.cxx: repeated
    /// multiplication, N >= 1).
    pub fn powered(&self, n: i32) -> GpMat {
        // OCCT gp_Mat::Powered: if (N == 1) return *this; if (N == 2)
        // return Multiplied(*this); if (N == 3) ... — generalized here as a
        // repeated product loop with identical term grouping per multiply.
        let mut result = *self;
        for _ in 1..n {
            result = result.multiplied_mat(self);
        }
        result
    }

    /// OCCT gp_Mat::Multiplied(const gp_Mat&) — matrix product (column
    /// convention: result(i,j) = sum_k this(i,k) * right(k,j)).
    pub fn multiplied_mat(&self, right: &GpMat) -> GpMat {
        let mut result = GpMat::identity();
        for (i, row) in result.mat.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = self.mat[i][0] * right.mat[0][j]
                    + self.mat[i][1] * right.mat[1][j]
                    + self.mat[i][2] * right.mat[2][j];
            }
        }
        result
    }

    /// OCCT gp_Mat::Multiplied(Scalar) — scale every coefficient.
    pub fn multiplied_scalar(&self, scalar: f64) -> GpMat {
        let mut result = *self;
        for row in result.mat.iter_mut() {
            for cell in row.iter_mut() {
                *cell *= scalar;
            }
        }
        result
    }

    /// OCCT gp_Mat::Added(const gp_Mat&) — coefficient sum.
    pub fn added(&self, right: &GpMat) -> GpMat {
        let mut result = *self;
        for (i, row) in result.mat.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell += right.mat[i][j];
            }
        }
        result
    }

    /// OCCT gp_Mat::Multiplied(const gp_XYZ&) — matrix x column-vector:
    /// result_i = sum_j this(i,j) * v_j (the `M * P` / `pnt * M` operator
    /// form used on gp_XYZ).
    pub fn multiplied_xyz(&self, the_xyz: DVec3) -> DVec3 {
        DVec3::new(
            self.mat[0][0] * the_xyz.x + self.mat[0][1] * the_xyz.y + self.mat[0][2] * the_xyz.z,
            self.mat[1][0] * the_xyz.x + self.mat[1][1] * the_xyz.y + self.mat[1][2] * the_xyz.z,
            self.mat[2][0] * the_xyz.x + self.mat[2][1] * the_xyz.y + self.mat[2][2] * the_xyz.z,
        )
    }

    /// OCCT gp_XYZ::Multiply(const gp_Mat&) — row-vector x matrix:
    /// result_i = sum_j v_j * this(j,i) (the in-place `pnt *= M` form).
    pub fn multiply_xyz_row(the_xyz: DVec3, m: &GpMat) -> DVec3 {
        DVec3::new(
            m.mat[0][0] * the_xyz.x + m.mat[1][0] * the_xyz.y + m.mat[2][0] * the_xyz.z,
            m.mat[0][1] * the_xyz.x + m.mat[1][1] * the_xyz.y + m.mat[2][1] * the_xyz.z,
            m.mat[0][2] * the_xyz.x + m.mat[1][2] * the_xyz.y + m.mat[2][2] * the_xyz.z,
        )
    }
}
