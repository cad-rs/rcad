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

use super::{IntVec, MatD, VecD};

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

/// OCCT SVD_Decompose file-local PYTHAG helper (math_Recipes.cxx L24-41,
/// anonymous namespace) — sqrt(a^2 + b^2) without destructive over/underflow.
fn pythag(a: f64, b: f64) -> f64 {
    let at = a.abs();
    let bt = b.abs();
    let mut ct = 0.0;
    if at > bt {
        ct = bt / at;
        ct = at * (1.0 + ct * ct).sqrt();
    } else if bt != 0.0 {
        ct = at / bt;
        ct = bt * (1.0 + ct * ct).sqrt();
    }
    ct
}

/// OCCT SVD_Decompose file-local SIGN macro (math_Recipes.cxx L48):
/// `((b) >= 0.0 ? fabs(a) : -fabs(a))`.
#[inline]
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

/// OCCT SVD_Decompose (math_Recipes.cxx L384-389) — the two-argument form
/// allocating the rv1 scratch vector internally.
pub fn svd_decompose(a: &mut MatD, w: &mut VecD, v: &mut MatD) -> i32 {
    let mut rv1 = VecD::new(a.n_cols());
    svd_decompose_with_rv1(a, w, v, &mut rv1)
}

/// OCCT SVD_Decompose (math_Recipes.cxx L391-710). Given a matrix
/// a(1..m, 1..n), computes its singular value decomposition a = u . diag(w) .
/// vT; `a` is replaced by U, the singular values are returned in `w`, and the
/// right singular vectors in `v`. Returns MATH_STATUS_OK or
/// MATH_STATUS_NO_CONVERGENCE.
pub fn svd_decompose_with_rv1(a: &mut MatD, w: &mut VecD, v: &mut MatD, rv1: &mut VecD) -> i32 {
    let mut flag: i32;
    let mut f: f64;
    let mut h: f64;
    let mut s: f64;
    let mut x: f64;
    let mut y: f64;
    let mut z: f64;
    let mut anorm = 0.0f64;
    let mut g = 0.0f64;
    let mut scale = 0.0f64;
    let m = a.n_rows() as i32;
    let n = a.n_cols() as i32;
    let mut l: i32 = 0;
    let mut nm: i32 = 0;

    for i in 1..=n {
        l = i + 1;
        rv1.set(i as usize, scale * g);
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if i <= m {
            for k in i..=m {
                let aki = a.get(k as usize, i as usize);
                if aki > 0.0 {
                    scale += aki;
                } else {
                    scale -= aki;
                }
            }
            if scale != 0.0 {
                for k in i..=m {
                    let v_ = a.get(k as usize, i as usize) / scale;
                    a.set(k as usize, i as usize, v_);
                    s += v_ * v_;
                }
                f = a.get(i as usize, i as usize);
                g = -sign(s.sqrt(), f);
                h = f * g - s;
                a.set(i as usize, i as usize, f - g);
                if i != n {
                    for j in l..=n {
                        s = 0.0;
                        for k in i..=m {
                            s += a.get(k as usize, i as usize) * a.get(k as usize, j as usize);
                        }
                        f = s / h;
                        for k in i..=m {
                            let v_ =
                                a.get(k as usize, j as usize) + f * a.get(k as usize, i as usize);
                            a.set(k as usize, j as usize, v_);
                        }
                    }
                }
                for k in i..=m {
                    let v_ = a.get(k as usize, i as usize) * scale;
                    a.set(k as usize, i as usize, v_);
                }
            }
        }
        w.set(i as usize, scale * g);
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if i <= m && i != n {
            for k in l..=n {
                let aik = a.get(i as usize, k as usize);
                if aik > 0.0 {
                    scale += aik;
                } else {
                    scale -= aik;
                }
            }
            if scale != 0.0 {
                for k in l..=n {
                    let v_ = a.get(i as usize, k as usize) / scale;
                    a.set(i as usize, k as usize, v_);
                    s += v_ * v_;
                }
                f = a.get(i as usize, l as usize);
                g = -sign(s.sqrt(), f);
                h = f * g - s;
                a.set(i as usize, l as usize, f - g);
                for k in l..=n {
                    let v_ = a.get(i as usize, k as usize) / h;
                    rv1.set(k as usize, v_);
                }
                if i != m {
                    for j in l..=m {
                        s = 0.0;
                        for k in l..=n {
                            s += a.get(j as usize, k as usize) * a.get(i as usize, k as usize);
                        }
                        for k in l..=n {
                            let v_ = a.get(j as usize, k as usize) + s * rv1.get(k as usize);
                            a.set(j as usize, k as usize, v_);
                        }
                    }
                }
                for k in l..=n {
                    let v_ = a.get(i as usize, k as usize) * scale;
                    a.set(i as usize, k as usize, v_);
                }
            }
        }
        let mut aw = w.get(i as usize);
        if aw < 0.0 {
            aw = -aw;
        }
        let mut ar = rv1.get(i as usize);
        if ar > 0.0 {
            ar = aw + ar;
        } else {
            ar = aw - ar;
        }
        if anorm < ar {
            anorm = ar;
        }
    }
    let mut i = n;
    while i >= 1 {
        if i < n {
            if g != 0.0 {
                for j in l..=n {
                    let v_ = (a.get(i as usize, j as usize) / a.get(i as usize, l as usize)) / g;
                    v.set(j as usize, i as usize, v_);
                }
                for j in l..=n {
                    s = 0.0;
                    for k in l..=n {
                        s += a.get(i as usize, k as usize) * v.get(k as usize, j as usize);
                    }
                    for k in l..=n {
                        let v_ =
                            v.get(k as usize, j as usize) + s * v.get(k as usize, i as usize);
                        v.set(k as usize, j as usize, v_);
                    }
                }
            }
            for j in l..=n {
                v.set(j as usize, i as usize, 0.0);
                v.set(i as usize, j as usize, 0.0);
            }
        }
        v.set(i as usize, i as usize, 1.0);
        g = rv1.get(i as usize);
        l = i;
        i -= 1;
    }
    let mut i = n;
    while i >= 1 {
        l = i + 1;
        g = w.get(i as usize);
        if i < n {
            for j in l..=n {
                a.set(i as usize, j as usize, 0.0);
            }
        }
        if g != 0.0 {
            g = 1.0 / g;
            if i != n {
                for j in l..=n {
                    s = 0.0;
                    for k in l..=m {
                        s += a.get(k as usize, i as usize) * a.get(k as usize, j as usize);
                    }
                    f = (s / a.get(i as usize, i as usize)) * g;
                    for k in i..=m {
                        let v_ = a.get(k as usize, j as usize) + f * a.get(k as usize, i as usize);
                        a.set(k as usize, j as usize, v_);
                    }
                }
            }
            for j in i..=m {
                let v_ = a.get(j as usize, i as usize) * g;
                a.set(j as usize, i as usize, v_);
            }
        } else {
            for j in i..=m {
                a.set(j as usize, i as usize, 0.0);
            }
        }
        let v_ = a.get(i as usize, i as usize) + 1.0; // OCCT L593: ++a(i, i);
        a.set(i as usize, i as usize, v_);
        i -= 1;
    }
    let mut k = n;
    while k >= 1 {
        let mut its: i32 = 1;
        while its <= 30 {
            flag = 1;
            l = k;
            while l >= 1 {
                nm = l - 1;
                if rv1.get(l as usize).abs() + anorm == anorm {
                    flag = 0;
                    break;
                }
                if w.get(nm as usize).abs() + anorm == anorm {
                    break;
                }
                l -= 1;
            }
            if flag != 0 {
                s = 1.0;
                for i in l..=k {
                    f = s * rv1.get(i as usize);
                    if f.abs() + anorm != anorm {
                        g = w.get(i as usize);
                        h = pythag(f, g);
                        w.set(i as usize, h);
                        h = 1.0 / h;
                        let c = g * h;
                        s = -f * h;
                        for j in 1..=m {
                            y = a.get(j as usize, nm as usize);
                            z = a.get(j as usize, i as usize);
                            a.set(j as usize, nm as usize, y * c + z * s);
                            a.set(j as usize, i as usize, z * c - y * s);
                        }
                    }
                }
            }
            z = w.get(k as usize);
            if l == k {
                if z < 0.0 {
                    w.set(k as usize, -z);
                    for j in 1..=n {
                        let v_ = -v.get(j as usize, k as usize);
                        v.set(j as usize, k as usize, v_);
                    }
                }
                break;
            }
            if its == 30 {
                return MATH_STATUS_NO_CONVERGENCE;
            }
            x = w.get(l as usize);
            nm = k - 1;
            y = w.get(nm as usize);
            g = rv1.get(nm as usize);
            h = rv1.get(k as usize);
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;

            let mut c = 1.0f64;
            s = 1.0f64;
            for j in l..=nm {
                let i = j + 1;
                g = rv1.get(i as usize);
                y = w.get(i as usize);
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1.set(j as usize, z);
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for jj in 1..=n {
                    x = v.get(jj as usize, j as usize);
                    z = v.get(jj as usize, i as usize);
                    v.set(jj as usize, j as usize, x * c + z * s);
                    v.set(jj as usize, i as usize, z * c - x * s);
                }
                z = pythag(f, h);
                w.set(j as usize, z);
                if z != 0.0 {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for jj in 1..=m {
                    y = a.get(jj as usize, j as usize);
                    z = a.get(jj as usize, i as usize);
                    a.set(jj as usize, j as usize, y * c + z * s);
                    a.set(jj as usize, i as usize, z * c - y * s);
                }
            }
            rv1.set(l as usize, 0.0);
            rv1.set(k as usize, f);
            w.set(k as usize, x);
            its += 1;
        }
        k -= 1;
    }
    MATH_STATUS_OK
}

/// OCCT SVD_Solve (math_Recipes.cxx L712-748). Solves u . diag(w) . vT . x = b
/// for x, where (u, w, v) are as returned by SVD_Decompose.
pub fn svd_solve(u: &MatD, w: &VecD, v: &MatD, b: &VecD, x: &mut VecD) {
    let mut s: f64;
    let m = u.n_rows();
    let n = u.n_cols();
    let mut tmp = VecD::new(n);

    for j in 1..=n {
        s = 0.0;
        if w.get(j) != 0.0 {
            for i in 1..=m {
                s += u.get(i, j) * b.get(i);
            }
            s /= w.get(j);
        }
        tmp.set(j, s);
    }
    for j in 1..=n {
        s = 0.0;
        for jj in 1..=n {
            s += v.get(j, jj) * tmp.get(jj);
        }
        x.set(j, s);
    }
}
