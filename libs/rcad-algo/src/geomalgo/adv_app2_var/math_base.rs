//! OCCT AdvApp2Var_MathBase (AdvApp2Var_MathBase.cxx) - the subset consumed
//! by the landed ApproxF2var engines and by AdvApp2Var_Patch::AddConstraints.
//!
//! Encoding note: the OCCT engine functions use Fortran 1-based storage with
//! `pointer -= offset` adjustments; rcad keeps the OCCT index formulas on
//! 0-based `&mut [f64]` buffers by subtracting the OCCT virtual origin
//! (`*_offset`) inside the flat-index expression, annotated per access.

use super::data_mlgdrtl::MLGDRTL_ROOTAB;
use super::sys_base;

/// OCCT AdvApp2Var_MathBase::pow__di (AdvApp2Var_MathBase.cxx L9530-9553) -
/// integer power, `1/x^|n|` for negative n.
pub fn pow__di(x: &f64, n: &i32) -> f64 {
    let mut result = 1.0e0;
    let absolute = if *n > 0 { *n } else { -*n };
    for _ii in 0..absolute {
        result *= *x;
    }
    if *n < 0 {
        result = 1.0e0 / result;
    }
    result
}

/// OCCT mfac_ (AdvApp2Var_MathBase.cxx L347-370) - fills `f[1..=n]` with the
/// factorials; the rcad slice is 1-based emulated (index 0 unused).
pub fn mfac_(f: &mut [f64], n: &i32) {
    f[1] = 1.;
    for i__ in 2..=*n {
        f[i__ as usize] = i__ as f64 * f[(i__ - 1) as usize];
    }
}

/// OCCT AdvApp2Var_MathBase::mmdrc11_ (AdvApp2Var_MathBase.cxx L2831-3013) -
/// successive derivatives of CURBE at parameters -1 and 1 (Horner scheme
/// generalized per Knuth, TAOCP Vol. 2 p. 423-425).
///
/// `courbe` is dimensioned (ncoeff, ndimen) Fortran-wise, `points` is
/// (2, iordre+1, ndimen) and `mfactab` the factorial workspace (1-based,
/// index 0 unused).
pub fn mmdrc11_(
    iordre: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    courbe: &[f64],
    points: &mut [f64],
    mfactab: &mut [f64],
) {
    let courbe_dim1 = *ncoeff;
    // OCCT: courbe_offset = courbe_dim1.
    let courbe_off = courbe_dim1 as usize;
    let points_dim2 = *iordre + 1;
    // OCCT: points_offset = (points_dim2 << 1) + 1.
    let points_off = (points_dim2 << 1) as usize + 1;
    let p = |i: i32, j: i32, nd: i32| -> usize {
        // OCCT: points[((i + nd * points_dim2) << 1) + j] - points_offset.
        (((i + nd * points_dim2) << 1) + j) as usize - points_off
    };
    let c = |i: i32, nd: i32| -> usize {
        // OCCT: courbe[i + nd * courbe_dim1] - courbe_offset.
        (i + nd * courbe_dim1) as usize - courbe_off
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 2 {
        sys_base::mgenmsg_("MMDRC11");
    }

    if *iordre < 0 || *ncoeff < 1 {
        return;
    }

    // ------------------- Initialization of table POINTS ------------------
    let ndgcb = *ncoeff - 1;
    for nd in 1..=*ndimen {
        points[p(0, 1, nd)] = courbe[c(ndgcb, nd)];
        points[p(0, 2, nd)] = courbe[c(ndgcb, nd)];
    }

    for nd in 1..=*ndimen {
        for j in 1..=*iordre {
            points[p(j, 1, nd)] = 0.;
            points[p(j, 2, nd)] = 0.;
        }
    }

    //    Calculation with parameter -1 and 1.
    for nd in 1..=*ndimen {
        for ndeg in 1..=ndgcb {
            for i__ in (1..=*iordre).rev() {
                points[p(i__, 1, nd)] = -points[p(i__, 1, nd)] + points[p(i__ - 1, 1, nd)];
                points[p(i__, 2, nd)] += points[p(i__ - 1, 2, nd)];
            }
            points[p(0, 1, nd)] = -points[p(0, 1, nd)] + courbe[c(ndgcb - ndeg, nd)];
            points[p(0, 2, nd)] += courbe[c(ndgcb - ndeg, nd)];
        }
    }

    // --------------------- Multiplication by factorial(I) ----------------
    if *iordre > 1 {
        // OCCT: mfac_(&mfactab[1], iordre).
        mfac_(mfactab, iordre);
        for nd in 1..=*ndimen {
            for i__ in 2..=*iordre {
                points[p(i__, 1, nd)] *= mfactab[i__ as usize];
                points[p(i__, 2, nd)] *= mfactab[i__ as usize];
            }
        }
    }

    if ibb >= 2 {
        sys_base::mgsomsg_("MMDRC11");
    }
}

/// OCCT AdvApp2Var_MathBase::mmrtptt_ (AdvApp2Var_MathBase.cxx L7921-8019) -
/// extracts the STRICTLY positive roots of the Legendre polynomial of degree
/// ndglgd (2 <= ndglgd <= 61) from the MLGDRTL block data (rootab).
///
/// `rtlegd` receives the RAW target address (the OCCT `&arr[k]` argument
/// maps to `&slice[(k - 1) as usize..]`); the OCCT body writes nsur2 doubles
/// at the received address (`&rtlegd[1]` after its `--rtlegd` adjustment
/// resolves to the raw pointer itself).
pub fn mmrtptt_(ndglgd: &mut i32, rtlegd: &mut [f64]) {
    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMRTPTT");
    }
    if *ndglgd < 2 {
        return;
    }

    let nsur2 = (*ndglgd / 2) as usize;
    let nmod2 = (*ndglgd % 2) as usize;

    let ilong = (nsur2 << 3) as i32;
    let ideb = nsur2 * (nsur2 - 1) / 2 + 1;
    // OCCT: mcrfill_(&ilong, &rootab[ideb + nmod2 * 465 - 1], &rtlegd[1]);
    // the source index is the 1-based block-data index ideb + nmod2*465.
    let src = &MLGDRTL_ROOTAB[ideb + nmod2 * 465 - 1..];
    sys_base::mcrfill_(ilong, src, &mut rtlegd[..nsur2]);

    if ibb >= 3 {
        sys_base::mgsomsg_("MMRTPTT");
    }
}
