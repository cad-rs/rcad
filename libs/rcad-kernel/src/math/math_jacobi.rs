//! OCCT math_Jacobi (TKMath) — eigenvalues and eigenvectors of a symmetric
//! square matrix by the cyclic Jacobi method.
//!
//! 1:1 translation of math_Jacobi.cxx (whole file) + the `Jacobi` and
//! `EigenSort` routines of math_Recipes.cxx (L60-90, L91-200).  Matrices and
//! vectors are the 1-based MatD/VecD.

use super::{MatD, VecD};

/// OCCT math_Status_OK.
const STATUS_OK: i32 = 0;
/// OCCT math_Status_NoConvergence.
const STATUS_NO_CONVERGENCE: i32 = 4;

/// OCCT ROTATE macro (math_Recipes.cxx).
#[inline]
fn rotate(a: &mut MatD, i: usize, j: usize, k: usize, l: usize, s: f64, tau: f64) {
    let g = a.get(i, j);
    let h = a.get(k, l);
    a.set(i, j, g - s * (h + g * tau));
    a.set(k, l, h + s * (g - h * tau));
}

/// OCCT static EigenSort (math_Recipes.cxx L60-90) — descending order.
fn eigen_sort(d: &mut VecD, v: &mut MatD) {
    let n = d.len();
    for i in 1..n {
        let mut k = i;
        let mut p = d.get(i);
        for j in (i + 1)..=n {
            if d.get(j) >= p {
                k = j;
                p = d.get(j);
            }
        }
        if k != i {
            d.set(k, d.get(i));
            d.set(i, p);
            for j in 1..=n {
                p = v.get(j, i);
                v.set(j, i, v.get(j, k));
                v.set(j, k, p);
            }
        }
    }
}

/// OCCT static Jacobi (math_Recipes.cxx L91-200) — cyclic Jacobi rotations;
/// returns math_Status_OK or math_Status_NoConvergence.
fn jacobi(a: &mut MatD, d: &mut VecD, v: &mut MatD, nrot: &mut i32) -> i32 {
    let n = a.n_rows() as i32;
    let mut b = VecD::new(n as usize);
    let mut z = VecD::new(n as usize);

    for ip in 1..=n {
        for iq in 1..=n {
            v.set(ip as usize, iq as usize, 0.0);
        }
        v.set(ip as usize, ip as usize, 1.0);
    }
    for ip in 1..=n {
        b.set(ip as usize, a.get(ip as usize, ip as usize));
        d.set(ip as usize, a.get(ip as usize, ip as usize));
        z.set(ip as usize, 0.0);
    }
    *nrot = 0;
    for i in 1..=50 {
        let mut sm = 0.0f64;
        for ip in 1..n {
            for iq in (ip + 1)..=n {
                sm += a.get(ip as usize, iq as usize).abs();
            }
        }
        if sm == 0.0 {
            eigen_sort(d, v);
            return STATUS_OK;
        }
        let tresh = if i < 4 {
            0.2 * sm / (n * n) as f64
        } else {
            0.0
        };
        for ip in 1..n {
            for iq in (ip + 1)..=n {
                let ipu = ip as usize;
                let iqu = iq as usize;
                let g = 100.0 * a.get(ipu, iqu).abs();
                if i > 4
                    && d.get(ipu).abs() + g == d.get(ipu).abs()
                    && d.get(iqu).abs() + g == d.get(iqu).abs()
                {
                    a.set(ipu, iqu, 0.0);
                } else if a.get(ipu, iqu).abs() > tresh {
                    let h = d.get(iqu) - d.get(ipu);
                    let t = if h.abs() + g == h.abs() {
                        a.get(ipu, iqu) / h
                    } else {
                        let theta = 0.5 * h / a.get(ipu, iqu);
                        let mut t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                        if theta < 0.0 {
                            t = -t;
                        }
                        t
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    let tau = s / (1.0 + c);
                    let h = t * a.get(ipu, iqu);
                    z.set(ipu, z.get(ipu) - h);
                    z.set(iqu, z.get(iqu) + h);
                    d.set(ipu, d.get(ipu) - h);
                    d.set(iqu, d.get(iqu) + h);
                    a.set(ipu, iqu, 0.0);
                    for j in 1..ip {
                        rotate(a, j as usize, ipu, j as usize, iqu, s, tau);
                    }
                    for j in (ip + 1)..iq {
                        rotate(a, ipu, j as usize, j as usize, iqu, s, tau);
                    }
                    for j in (iq + 1)..=n {
                        rotate(a, ipu, j as usize, iqu, j as usize, s, tau);
                    }
                    for j in 1..=n {
                        rotate(v, j as usize, ipu, j as usize, iqu, s, tau);
                    }
                    *nrot += 1;
                }
            }
        }
        for ip in 1..=n {
            b.set(ip as usize, b.get(ip as usize) + z.get(ip as usize));
            d.set(ip as usize, b.get(ip as usize));
            z.set(ip as usize, 0.0);
        }
    }
    eigen_sort(d, v);
    STATUS_NO_CONVERGENCE
}

/// OCCT math_Jacobi — eigen decomposition of a symmetric square matrix.
#[derive(Debug, Clone)]
pub struct MathJacobi {
    done: bool,
    aa: MatD,
    nb_rotations: i32,
    eigen_values: VecD,
    eigen_vectors: MatD,
}

impl MathJacobi {
    /// OCCT math_Jacobi(A) (math_Jacobi.cxx L27-40).  A must be square;
    /// panics (math_NotSquare) otherwise.
    pub fn new(a: &MatD) -> Self {
        assert!(
            a.n_rows() == a.n_cols(),
            "math_NotSquare: math_Jacobi requires a square matrix"
        );
        let n = a.n_rows();
        let mut aa = MatD::new(n, n);
        for i in 1..=n {
            for j in 1..=n {
                aa.set(i, j, a.get(i, j));
            }
        }
        let mut eigen_values = VecD::new(n);
        let mut eigen_vectors = MatD::new(n, n);
        let mut nb_rotations = 0i32;
        let error = jacobi(
            &mut aa,
            &mut eigen_values,
            &mut eigen_vectors,
            &mut nb_rotations,
        );
        MathJacobi {
            done: error == 0,
            aa,
            nb_rotations,
            eigen_values,
            eigen_vectors,
        }
    }

    /// OCCT IsDone (hxx).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Values() — the eigenvalues, in descending order.
    pub fn values(&self) -> &VecD {
        &self.eigen_values
    }

    /// OCCT Vectors() — the eigenvectors (column i belongs to value i).
    pub fn vectors(&self) -> &MatD {
        &self.eigen_vectors
    }

    /// OCCT math_Jacobi::Value(Num) (math_Jacobi.lxx) — the eigenvalue
    /// number Num.  Eigenvalues are in the range (1..n).
    pub fn value(&self, num: usize) -> f64 {
        self.eigen_values.get(num)
    }

    /// OCCT math_Jacobi::Vector(Num, V) (math_Jacobi.lxx) — returns the
    /// eigenvector V of number Num (column Num of the eigenvector matrix).
    /// Eigenvectors are in the range (1..n).
    pub fn vector(&self, num: usize) -> VecD {
        let n = self.eigen_vectors.n_rows();
        let mut v = VecD::new(n);
        for j in 1..=n {
            v.set(j, self.eigen_vectors.get(j, num));
        }
        v
    }

    /// OCCT AA() — the input matrix (unchanged copy).
    pub fn aa(&self) -> &MatD {
        &self.aa
    }

    /// OCCT NbRotations() (hxx).
    pub fn nb_rotations(&self) -> i32 {
        self.nb_rotations
    }
}
