// OCCT math_Recipes (TKMath/math/math_Recipes.hxx/.cxx) — the numerical
// recipes used by the AppParCurves/AppDef approximation templates.
//
// Ported so far: DACTCL_Decompose + DACTCL_Solve (math_Recipes.cxx L750-901)
// — the symmetric banded LDLT decomposition / solve pair called by
// AppParCurves_LeastSquare.gxx (Perform / Error). LU_Decompose / LU_Solve
// already live in math_gauss.rs. The remaining recipes (SVD, Jacobi, FFT...)
// are added on demand by their consuming alignment units.
//
// OCCT status codes (math_Recipes.hxx L14-19).

/// OCCT math_Status_UserAborted.
pub const MATH_STATUS_USER_ABORTED: i32 = -1;
/// OCCT math_Status_OK.
pub const MATH_STATUS_OK: i32 = 0;
/// OCCT math_Status_SingularMatrix.
pub const MATH_STATUS_SINGULAR_MATRIX: i32 = 1;
/// OCCT math_Status_ArgumentError.
pub const MATH_STATUS_ARGUMENT_ERROR: i32 = 2;
/// OCCT math_Status_NoConvergence.
pub const MATH_STATUS_NO_CONVERGENCE: i32 = 3;

use super::{IntVec, VecD};

/// OCCT DACTCL_Decompose (math_Recipes.cxx L750-830). Given a SYMMETRIC
/// matrix `a`, computes its LU decomposition. `a` is given through a vector
/// of its non-zero components of the upper triangular matrix; `indx` is the
/// index vector of the diagonal elements of `a`; `a` is replaced by its LU
/// decomposition. The range of the matrix is n = indx.Length(), and
/// a.Length() = indx(n). OCCT MinPivot default: 1.0e-20.
pub fn dactcl_decompose(a: &mut VecD, indx: &IntVec, min_pivot: f64) -> i32 {
    let neq = indx.len() as i32;
    let mut jr: i32 = 0;
    let mut j: i32 = 1;
    while j <= neq {
        let mut diag = false;
        let jd = indx.get(j as usize);
        let jh = jd - jr;
        let is = j - jh + 2;
        if jh - 2 == 0 {
            diag = true;
        }
        if jh - 2 > 0 {
            let ie = j - 1;
            let mut k = jr + 2;
            let mut id = indx.get((is - 1) as usize);
            // Reduction des coefficients non diagonaux:
            // =========================================
            let mut i = is;
            while i <= ie {
                let ir = id;
                id = indx.get(i as usize);
                let mut ih = id - ir - 1;
                let mh = i - is + 1;
                if ih > mh {
                    ih = mh;
                }
                if ih > 0 {
                    let mut dot = 0.0;
                    let idot1 = k - ih - 1;
                    let idot2 = id - ih - 1;
                    for idot in 1..=ih {
                        dot += a.get((idot1 + idot) as usize) * a.get((idot2 + idot) as usize);
                    }
                    let v = a.get(k as usize) - dot;
                    a.set(k as usize, v);
                }
                k += 1;
                i += 1;
            }
            diag = true;
        }

        if diag {
            // Reduction des coefficients diagonaux:
            // =====================================
            let ir = jr + 1;
            let ie = jd - 1;
            let k = j - jd;
            let mut i = ir;
            while i <= ie {
                let id = indx.get((k + i) as usize);
                let mut aa = a.get(id as usize);
                if aa < 0.0 {
                    aa = -aa;
                }
                if aa <= min_pivot {
                    return MATH_STATUS_SINGULAR_MATRIX;
                }
                let d = a.get(i as usize);
                let v = d / a.get(id as usize);
                a.set(i as usize, v);
                let w = a.get(jd as usize) - d * a.get(i as usize);
                a.set(jd as usize, w);
                i += 1;
            }
        }
        jr = jd;
        j += 1;
    }
    MATH_STATUS_OK
}

/// OCCT DACTCL_Solve (math_Recipes.cxx L832-901). Solves a * x = b for a
/// vector x and a matrix `a` coming from DACTCL_Decompose; `indx` is the
/// same vector as in DACTCL_Decompose. The vector `b` is replaced by the
/// vector solution x. OCCT MinPivot default: 1.0e-20.
pub fn dactcl_solve(a: &VecD, b: &mut VecD, indx: &IntVec, min_pivot: f64) -> i32 {
    let neq = indx.len() as i32;
    let mut jr: i32 = 0;
    let mut j: i32 = 1;
    while j <= neq {
        let jd = indx.get(j as usize);
        let jh = jd - jr;
        let is = j - jh + 2;

        // Reduction du second membre:
        // ===========================
        let mut dot = 0.0;
        let idot1 = jr;
        let idot2 = is - 2;
        let jh1 = jh - 1;
        for idot in 1..=jh1 {
            dot += a.get((idot1 + idot) as usize) * b.get((idot2 + idot) as usize);
        }
        let v = b.get(j as usize) - dot;
        b.set(j as usize, v);

        jr = jd;
        j += 1;
    }

    // Division par les pivots diagonaux:
    // ==================================
    let mut i: i32 = 1;
    while i <= neq {
        let id = indx.get(i as usize);
        let mut aa = a.get(id as usize);
        if aa < 0.0 {
            aa = -aa;
        }
        if aa <= min_pivot {
            return MATH_STATUS_SINGULAR_MATRIX;
        }
        let v = b.get(i as usize) / a.get(id as usize);
        b.set(i as usize, v);
        i += 1;
    }

    // Substitution arriere:
    // =====================
    let mut jd = indx.get(neq as usize);
    let mut j = neq - 1;
    while j > 0 {
        let d = b.get((j + 1) as usize);
        let jr = indx.get(j as usize);
        if jd - jr > 1 {
            let is = j - jd + jr + 2;
            let k = jr - is + 1;
            for i in is..=j {
                let v = b.get(i as usize) - a.get((i + k) as usize) * d;
                b.set(i as usize, v);
            }
        }
        jd = jr;
        j -= 1;
    }
    MATH_STATUS_OK
}
