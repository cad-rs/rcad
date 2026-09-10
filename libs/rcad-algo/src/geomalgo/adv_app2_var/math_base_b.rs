//! OCCT AdvApp2Var_MathBase (AdvApp2Var_MathBase.cxx) - part B: the MathBase
//! leaves consumed by the ApproxF2var engines (mzsnorm_, mmveps3_, mmeps1_,
//! mmapcmp_, mmaper0/2/4/6/x, mmfmca8_, mmfmca9_, mmfmtb1_, mmjacan_,
//! mmjaccv_, mmmpocur_, mmtrpjj_ + mmtrpj0/2/4/6) and the MPRCSN common
//! writer mmwprcs_.
//!
//! Encoding note: as in part A ([`super::math_base`]), the OCCT Fortran
//! 1-based index formulas are kept on 0-based slices by subtracting the OCCT
//! virtual origin inside the flat-index expression (annotated per function
//! as `*_offset`).  `&arr[k]` (1-based k) passes as `&arr[(k - 1) as
//! usize..]`.

use std::sync::atomic::AtomicI32;
use std::sync::atomic::AtomicU64;
use std::sync::atomic::Ordering;

use super::approxf2var_b::XMAX2;
use super::approxf2var_b::XMAX4;
use super::approxf2var_b::XMAX6;
use super::data_mmjcobi::PLGCAN;
use super::sys_base;

// ---------------------------------------------------------------------------
// COMMON MPRCSN (AdvApp2Var_MathBase.cxx L219-226) and its writer mmwprcs_.
//
// OCCT keeps `static struct { double eps1, eps2, eps3, eps4; int niterm,
// niterr; } mmprcsn_;` as a translation-unit static initialized once before
// first public use via mmwprcs_ (initMathBasePrecision, L24-58: EPS1=1e-9,
// EPS2=1e-8, EPS3=1e-9, EPS4=1e-4, ITERM=8, ITERR=40).  rcad encodes the
// struct in atomics (f64 as bit patterns) with the same one-time
// initialization semantics.
// ---------------------------------------------------------------------------

static MPRCSN_EPS1: AtomicU64 = AtomicU64::new(0);
static MPRCSN_EPS2: AtomicU64 = AtomicU64::new(0);
static MPRCSN_EPS3: AtomicU64 = AtomicU64::new(0);
static MPRCSN_EPS4: AtomicU64 = AtomicU64::new(0);
static MPRCSN_NITERM: AtomicI32 = AtomicI32::new(0);
static MPRCSN_NITERR: AtomicI32 = AtomicI32::new(0);
static MPRCSN_INIT: AtomicU64 = AtomicU64::new(0);

/// One-time static initialization trigger - OCCT `THE_STMAT_LIB_INIT`
/// (AdvApp2Var_MathBase.cxx L51-56) / initMathBasePrecision (L24-47).
fn init_math_base_precision() {
    if MPRCSN_INIT.load(Ordering::Acquire) == 0 {
        // constexpr double THE_EPSILON_1 = 1.0e-9; ... (cxx L25-30).
        mmwprcs_(
            &mut 1.0e-9,
            &mut 1.0e-8,
            &mut 1.0e-9,
            &mut 1.0e-4,
            &mut 8,
            &mut 40,
        );
        MPRCSN_INIT.store(1, Ordering::Release);
    }
}

fn mprcsn_eps1() -> f64 {
    init_math_base_precision();
    f64::from_bits(MPRCSN_EPS1.load(Ordering::Acquire))
}

fn mprcsn_eps3() -> f64 {
    init_math_base_precision();
    f64::from_bits(MPRCSN_EPS3.load(Ordering::Acquire))
}

/// OCCT AdvApp2Var_MathBase::mmwprcs_ (AdvApp2Var_MathBase.cxx L9427-9526) -
/// write access to COMMON MPRCSN.
pub fn mmwprcs_(
    epsil1: &mut f64,
    epsil2: &mut f64,
    epsil3: &mut f64,
    epsil4: &mut f64,
    niter1: &mut i32,
    niter2: &mut i32,
) {
    // OCCT: mmprcsn_.eps1 = *epsil1; ... (cxx L9520-9525).
    MPRCSN_EPS1.store(epsil1.to_bits(), Ordering::Release);
    MPRCSN_EPS2.store(epsil2.to_bits(), Ordering::Release);
    MPRCSN_EPS3.store(epsil3.to_bits(), Ordering::Release);
    MPRCSN_EPS4.store(epsil4.to_bits(), Ordering::Release);
    MPRCSN_NITERM.store(*niter1, Ordering::Release);
    MPRCSN_NITERR.store(*niter2, Ordering::Release);
}

/// OCCT AdvApp2Var_MathBase::mmeps1_ (AdvApp2Var_MathBase.cxx L3385-3478) -
/// extraction of EPS1 from COMMON MPRCSN (spatial zero, 1.0e-9).
pub fn mmeps1_(epsilo: &mut f64) {
    // OCCT: *epsilo = mmprcsn_.eps1 (cxx L3463).
    *epsilo = mprcsn_eps1();
}

/// OCCT AdvApp2Var_MathBase::mzsnorm_ (AdvApp2Var_MathBase.cxx L10503-10619) -
/// euclidean norm of a vector; the term of strongest absolute value is
/// factorized to limit overflow risks.  `vecteu` is 1-based emulated (the
/// OCCT body adjusts `--vecteu`).
pub fn mzsnorm_(ndimen: &mut i32, vecteu: &[f64]) -> f64 {
    // OCCT: --vecteu (1-based emulated slice, index 0 unused).
    let v = |i: i32| -> usize {
        // OCCT: vecteu[i] (1-based).
        (i - 1) as usize
    };

    // ___ Find the strongest absolute value term
    let mut irmax = 1i32;
    for i__ in 2..=*ndimen {
        if vecteu[v(irmax)].abs() < vecteu[v(i__)].abs() {
            irmax = i__;
        }
    }

    // ___ Calculate the norme
    let ret_val;
    if vecteu[v(irmax)].abs() < 1. {
        let mut xsom = 0.;
        for i__ in 1..=*ndimen {
            xsom += vecteu[v(i__)] * vecteu[v(i__)];
        }
        ret_val = xsom.sqrt();
    } else {
        let mut xsom = 0.;
        for i__ in 1..=*ndimen {
            if i__ == irmax {
                xsom += 1.;
            } else {
                let d__1 = vecteu[v(i__)] / vecteu[v(irmax)];
                xsom += d__1 * d__1;
            }
        }
        ret_val = vecteu[v(irmax)].abs() * xsom.sqrt();
    }
    ret_val
}

/// OCCT AdvApp2Var_MathBase::mmveps3_ (AdvApp2Var_MathBase.cxx L9171-9261) -
/// extraction of EPS3 from COMMON MPRCSN (10**-9).
pub fn mmveps3_(eps03: &mut f64) {
    let ibb = sys_base::mnfndeb_();
    if ibb >= 5 {
        sys_base::mgenmsg_("MMEPS1  ");
    }
    // OCCT: *eps03 = mmprcsn_.eps3 (cxx L9256).
    *eps03 = mprcsn_eps3();
}

/// OCCT AdvApp2Var_MathBase::mmapcmp_ (AdvApp2Var_MathBase.cxx L374-486) -
/// compression of CRVOLD into the even terms CRVNEW(*,0,*) and uneven terms
/// CRVNEW(*,1,*).
pub fn mmapcmp_(ndim: &mut i32, ncofmx: &mut i32, ncoeff: &mut i32, crvold: &[f64], crvnew: &mut [f64]) {
    // OCCT: crvold_dim1 = *ncofmx; crvold_offset = crvold_dim1.
    let crvold_dim1 = *ncofmx;
    let crvold_off = crvold_dim1;
    // OCCT: crvnew_dim1 = (*ncoeff - 1) / 2 + 1; crvnew_offset = dim1 << 1.
    let crvnew_dim1 = ((*ncoeff - 1) / 2 + 1) as usize;
    let crvnew_off = crvnew_dim1 << 1;

    let co = |i: i32, nd: i32| -> usize {
        // OCCT: crvold[i + nd * crvold_dim1] - crvold_offset.
        ((i + nd * crvold_dim1) - crvold_off) as usize
    };
    let cn = |idg: i32, col: i32| -> usize {
        // OCCT: crvnew[idg + col * crvnew_dim1] - crvnew_offset.
        ((idg + col * crvnew_dim1 as i32) - crvnew_off as i32) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMAPCMP");
    }

    let ndegre = *ncoeff - 1;
    for nd in 1..=*ndim {
        let mut ipair = 0i32;
        for idg in 0..=ndegre / 2 {
            crvnew[cn(idg, nd << 1)] = crvold[co(ipair, nd)];
            ipair += 2;
        }
        if ndegre < 1 {
            // OCCT goto L400 (next nd).
            continue;
        }
        let mut impair = 1i32;
        for idg in 0..=(ndegre - 1) / 2 {
            crvnew[cn(idg, (nd << 1) + 1)] = crvold[co(impair, nd)];
            impair += 2;
        }
    }

    if ibb >= 3 {
        sys_base::mgsomsg_("MMAPCMP");
    }
}

/// OCCT mmaper0_ (AdvApp2Var_MathBase.cxx L487-594) - max error when only
/// the first NCFNEW coefficients of a Legendre-basis curve are preserved.
pub fn mmaper0_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    crvlgd: &[f64],
    ncfnew: &mut i32,
    ycvmax: &mut [f64],
    errmax: &mut f64,
) {
    // OCCT: crvlgd_dim1 = *ncofmx; crvlgd_offset = crvlgd_dim1 + 1.
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        // OCCT: crvlgd[i + nd * crvlgd_dim1] - crvlgd_offset.
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    for ii in 1..=*ndimen {
        // OCCT: ycvmax[ii] (1-based).
        ycvmax[(ii - 1) as usize] = 0.;
    }

    // ------ Minimum that can be reached : Stop at 1 or NCFNEW ------
    let mut ncut = 1;
    if *ncfnew >= ncut {
        ncut = *ncfnew + 1;
    }

    // ------ Elimination of high degree coefficients: NCUT --> NCOEFF ------
    for ii in ncut..=*ncoeff {
        //   Factor of renormalization (Maximum of Li(t)).
        let mut bidon = ((ii - 1) as f64 * 2. + 1.) / 2.;
        bidon = bidon.sqrt();
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvlgd[cl(ii, nd)].abs() * bidon;
        }
    }

    // ------ The error is the norm of the vector error ------
    *errmax = mzsnorm_(ndimen, ycvmax);
}

/// OCCT mmaper2_ (AdvApp2Var_MathBase.cxx L595-737) - same as mmaper0_ in
/// the Jacobi basis of order 2.  OCCT carries an inline `xmaxj[57]` static
/// identical to mma2jmx_'s XMAX2 (approxf2var_b.rs).
pub fn mmaper2_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    crvjac: &[f64],
    ncfnew: &mut i32,
    ycvmax: &mut [f64],
    errmax: &mut f64,
) {
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    for ii in 1..=*ndimen {
        ycvmax[(ii - 1) as usize] = 0.;
    }

    // ------ Min. Degree that can be attained : Stop at 3 or NCFNEW ------
    let idec = 3;
    let ncut = idec.max(*ncfnew + 1);

    for ii in ncut..=*ncoeff {
        //   Factor of renormalization (maximum of (1-t^2)*Ji(t)).
        let bidon = XMAX2[(ii - idec) as usize];
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvjac[cl(ii, nd)].abs() * bidon;
        }
    }

    *errmax = mzsnorm_(ndimen, ycvmax);
}

/// OCCT mmaper4_ (AdvApp2Var_MathBase.cxx L738-873) - Jacobi basis of
/// order 4.  OCCT carries an inline `xmaxj[55]` static identical to
/// mma2jmx_'s XMAX4 (approxf2var_b.rs).
pub fn mmaper4_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    crvjac: &[f64],
    ncfnew: &mut i32,
    ycvmax: &mut [f64],
    errmax: &mut f64,
) {
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    for ii in 1..=*ndimen {
        ycvmax[(ii - 1) as usize] = 0.;
    }

    // ------ Min. Degree that can be attained : Stop at 5 or NCFNEW ------
    let idec = 5;
    let ncut = idec.max(*ncfnew + 1);

    for ii in ncut..=*ncoeff {
        let bidon = XMAX4[(ii - idec) as usize];
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvjac[cl(ii, nd)].abs() * bidon;
        }
    }

    *errmax = mzsnorm_(ndimen, ycvmax);
}

/// OCCT mmaper6_ (AdvApp2Var_MathBase.cxx L874-1008) - Jacobi basis of
/// order 6.  OCCT carries an inline `xmaxj[53]` static identical to
/// mma2jmx_'s XMAX6 (approxf2var_b.rs).
pub fn mmaper6_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    crvjac: &[f64],
    ncfnew: &mut i32,
    ycvmax: &mut [f64],
    errmax: &mut f64,
) {
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    for ii in 1..=*ndimen {
        ycvmax[(ii - 1) as usize] = 0.;
    }

    // ------ Min Degree that can be attained : Stop at 7 or NCFNEW ------
    let idec = 7;
    let ncut = idec.max(*ncfnew + 1);

    for ii in ncut..=*ncoeff {
        let bidon = XMAX6[(ii - idec) as usize];
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvjac[cl(ii, nd)].abs() * bidon;
        }
    }

    *errmax = mzsnorm_(ndimen, ycvmax);
}

/// OCCT AdvApp2Var_MathBase::mmaperx_ (AdvApp2Var_MathBase.cxx L1009-1062) -
/// max error when only the first NCFNEW coefficients of a Jacobi-basis curve
/// of order IORDRE are preserved; dispatches on jord = 2*(IORDRE+1).
pub fn mmaperx_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    iordre: &mut i32,
    crvjac: &[f64],
    ncfnew: &mut i32,
    ycvmax: &mut [f64],
    errmax: &mut f64,
    iercod: &mut i32,
) {
    *iercod = 0;
    // --> Order of Jacobi polynoms
    let jord = (*iordre + 1) << 1;

    // OCCT: the crvjac sub-pointer &crvjac[crvjac_offset] resolves to the
    // raw address (the rcad slice is the physical buffer).
    if jord == 0 {
        mmaper0_(ncofmx, ndimen, ncoeff, crvjac, ncfnew, ycvmax, errmax);
    } else if jord == 2 {
        mmaper2_(ncofmx, ndimen, ncoeff, crvjac, ncfnew, ycvmax, errmax);
    } else if jord == 4 {
        mmaper4_(ncofmx, ndimen, ncoeff, crvjac, ncfnew, ycvmax, errmax);
    } else if jord == 6 {
        mmaper6_(ncofmx, ndimen, ncoeff, crvjac, ncfnew, ycvmax, errmax);
    } else {
        *iercod = 1;
    }
}

/// OCCT AdvApp2Var_MathBase::mmfmca8_ (AdvApp2Var_MathBase.cxx L3745-3871) -
/// expansion of a table containing only the most important things into a
/// greater data table (decompression NDIMAX<>NDIMEN / =NDIMEN cases).
pub fn mmfmca8_(
    ndimen: &i32,
    ncoefu: &i32,
    ncoefv: &i32,
    ndimax: &i32,
    ncfumx: &i32,
    ncfvmx: &i32,
    tabini: &[f64],
    tabres: &mut [f64],
) {
    let _ = ncfvmx;
    // OCCT: tabini_dim1 = *ndimen; tabini_dim2 = *ncoefu;
    //       tabini_offset = tabini_dim1 * (tabini_dim2 + 1) + 1.
    let tabini_dim1 = *ndimen;
    let tabini_dim2 = *ncoefu;
    let tabini_off = tabini_dim1 * (tabini_dim2 + 1) + 1;
    // OCCT: tabres_dim1 = *ndimax; tabres_dim2 = *ncfumx;
    //       tabres_offset = tabres_dim1 * (tabres_dim2 + 1) + 1.
    let tabres_dim1 = *ndimax;
    let tabres_dim2 = *ncfumx;
    let tabres_off = tabres_dim1 * (tabres_dim2 + 1) + 1;

    let ti = |i: i32, j: i32, k: i32| -> usize {
        // OCCT: tabini[i + (j + k * tabini_dim2) * tabini_dim1] - offset.
        ((i + (j + k * tabini_dim2) * tabini_dim1) - tabini_off) as usize
    };
    let tr = |i: i32, j: i32, k: i32| -> usize {
        // OCCT: tabres[i + (j + k * tabres_dim2) * tabres_dim1] - offset.
        ((i + (j + k * tabres_dim2) * tabres_dim1) - tabres_off) as usize
    };

    if *ndimax != *ndimen {
        // ---------------- decompression NDIMAX<>NDIMEN ----------------
        for k in (1..=*ncoefv).rev() {
            for j in (1..=*ncoefu).rev() {
                for i__ in (1..=*ndimen).rev() {
                    tabres[tr(i__, j, k)] = tabini[ti(i__, j, k)];
                }
            }
        }
        return; // OCCT goto L9999.
    }

    // ---------------- decompression NDIMAX=NDIMEN ----------------
    if *ncoefu != *ncfumx {
        let ilong = (*ndimen << 3) * *ncoefu;
        for k in (1..=*ncoefv).rev() {
            // OCCT: &tabini[(k * tabini_dim2 + 1) * tabini_dim1 + 1]
            //   resolves to raw index (k * tabini_dim2 + 1) * tabini_dim1
            //   + 1 - tabini_offset = k * tabini_dim2 * tabini_dim1.
            let src = (k * tabini_dim2 * tabini_dim1) as usize;
            let dst = (k * tabres_dim2 * tabres_dim1) as usize;
            sys_base::mcrfill_(ilong, &tabini[src..], &mut tabres[dst..]);
        }
        return; // OCCT goto L9999.
    }

    // ---------------- decompression NDIMAX=NDIMEN,NCOEFU=NCFUMX ----------
    let ilong = (*ndimen << 3) * *ncoefu * *ncoefv;
    sys_base::mcrfill_(ilong, tabini, tabres);
}

/// OCCT AdvApp2Var_MathBase::mmfmca9_ (AdvApp2Var_MathBase.cxx L3873-4008) -
/// compression of a data table into a table containing only the main data
/// (the input table is not removed).
pub fn mmfmca9_(
    ndimax: &i32,
    ncfumx: &i32,
    ncfvmx: &i32,
    ndimen: &i32,
    ncoefu: &i32,
    ncoefv: &i32,
    tabini: &[f64],
    tabres: &mut [f64],
) {
    let _ = ncfvmx;
    // OCCT: tabini_dim1 = *ndimax; tabini_dim2 = *ncfumx;
    //       tabini_offset = tabini_dim1 * (tabini_dim2 + 1) + 1.
    let tabini_dim1 = *ndimax;
    let tabini_dim2 = *ncfumx;
    let tabini_off = tabini_dim1 * (tabini_dim2 + 1) + 1;
    // OCCT: tabres_dim1 = *ndimen; tabres_dim2 = *ncoefu;
    //       tabres_offset = tabres_dim1 * (tabres_dim2 + 1) + 1.
    let tabres_dim1 = *ndimen;
    let tabres_dim2 = *ncoefu;
    let tabres_off = tabres_dim1 * (tabres_dim2 + 1) + 1;

    let ti = |i: i32, j: i32, k: i32| -> usize {
        ((i + (j + k * tabini_dim2) * tabini_dim1) - tabini_off) as usize
    };
    let tr = |i: i32, j: i32, k: i32| -> usize {
        ((i + (j + k * tabres_dim2) * tabres_dim1) - tabres_off) as usize
    };

    if *ndimen != *ndimax {
        // ---------------- compression NDIMEN<>NDIMAX ----------------
        for k in 1..=*ncoefv {
            for j in 1..=*ncoefu {
                for i__ in 1..=*ndimen {
                    tabres[tr(i__, j, k)] = tabini[ti(i__, j, k)];
                }
            }
        }
        return; // OCCT goto L9999.
    }

    // ---------------- compression NDIMEN=NDIMAX ----------------
    if *ncoefu != *ncfumx {
        let ilong = (*ndimen << 3) * *ncoefu;
        for k in 1..=*ncoefv {
            let src = (k * tabini_dim2 * tabini_dim1) as usize;
            let dst = (k * tabres_dim2 * tabres_dim1) as usize;
            sys_base::mcrfill_(ilong, &tabini[src..], &mut tabres[dst..]);
        }
        return; // OCCT goto L9999.
    }

    // ---------------- compression NDIMEN=NDIMAX,NCOEFU=NCFUMX -----------
    let ilong = (*ndimen << 3) * *ncoefu * *ncoefv;
    sys_base::mcrfill_(ilong, tabini, tabres);
}

/// OCCT AdvApp2Var_MathBase::mmfmtb1_ (AdvApp2Var_MathBase.cxx L4412-4553) -
/// transposition of a rectangular table (T1(i,j) loaded in T2(j,i)).
pub fn mmfmtb1_(
    maxsz1: &mut i32,
    table1: &[f64],
    isize1: &mut i32,
    jsize1: &mut i32,
    maxsz2: &mut i32,
    table2: &mut [f64],
    isize2: &mut i32,
    jsize2: &mut i32,
    iercod: &mut i32,
) {
    let c__8: i32 = 8;

    // OCCT: table1_dim1 = *maxsz1; table1_offset = table1_dim1 + 1.
    let table1_dim1 = *maxsz1;
    let table1_off = table1_dim1 + 1;
    // OCCT: table2_dim1 = *maxsz2; table2_offset = table2_dim1 + 1.
    let table2_dim1 = *maxsz2;
    let table2_off = table2_dim1 + 1;

    let t1 = |i: i32, j: i32| -> usize {
        // OCCT: table1[i + j * table1_dim1] - table1_offset.
        ((i + j * table1_dim1) - table1_off) as usize
    };

    *iercod = 0;
    if *isize1 > *maxsz1 || *jsize1 > *maxsz2 {
        // OCCT goto L9100.
        *iercod = 1;
        finish_mmfmtb1_(iercod);
        return;
    }

    let mut work: Vec<f64> = Vec::new();
    let mut iofst: isize = 0;
    let mut isize_ = *maxsz2 * *isize1;
    let mut ier = 0;
    let mut iunit = c__8;
    sys_base::mcrrqst_(&mut iunit, &mut isize_, &mut work, &mut iofst, &mut ier);
    if ier > 0 {
        // OCCT goto L9200.
        *iercod = 2;
        finish_mmfmtb1_(iercod);
        return;
    }

    //         DO NOT BE AFRAID OF CRUSHING.
    for ii in 1..=*isize1 {
        let iipt = ((ii - 1) * *maxsz2) as isize + iofst;
        for jj in 1..=*jsize1 {
            let jjpt = iipt + (jj - 1) as isize;
            work[jjpt as usize] = table1[t1(ii, jj)];
        }
    }
    let ilong = isize_ << 3;
    // OCCT: &work[iofst] resolves to the raw workspace address (index 0);
    // &table2[table2_offset] resolves to the raw table address (index 0).
    sys_base::mcrfill_(ilong, &work[..], table2);

    // -------------- The number of elements of TABLE2 is returned ----------
    let ii = *isize1;
    *isize2 = *jsize1;
    *jsize2 = ii;

    finish_mmfmtb1_(iercod);
}

/// OCCT mmfmtb1_ L9999 tail (mcrdelt_ + error overwrite).
fn finish_mmfmtb1_(iercod: &mut i32) {
    let c__8: i32 = 8;
    // OCCT keeps `ier` at 0 across the L9999 fall-through; mcrdelt_ never
    // fails in the rcad encoding, matching OCCT's ier > 0 guard.
    let mut ier = 0;
    let mut isize_ = 0;
    let mut iofst: isize = 0;
    let mut t: Vec<f64> = Vec::new();
    let mut iunit = c__8;
    sys_base::mcrdelt_(&mut iunit, &mut isize_, &mut t, &mut iofst, &mut ier);
    if ier > 0 {
        *iercod = 2;
    }
}

/// OCCT AdvApp2Var_MathBase::mmjacan_ (AdvApp2Var_MathBase.cxx L5660-5791) -
/// transfer Jacobi normalized to canonic [-1,1]; tables ranked by even then
/// uneven degree.  `poljac` is a 0-based raw buffer (dim ndeg/2+1 per half),
/// `polcan` is the 0-based raw sub-buffer handed by the caller.
pub fn mmjacan_(ideriv: &i32, ndeg: &mut i32, poljac: &[f64], polcan: &mut [f64]) {
    // OCCT: poljac_dim1 = *ndeg / 2 + 1.
    let poljac_dim1 = *ndeg / 2 + 1;

    let ibb = sys_base::mnfndeb_();
    if ibb >= 5 {
        sys_base::mgenmsg_("MMJACAN");
    }

    // ----------------- Expression of terms of even degree ----------------
    for i__ in 0..=*ndeg / 2 {
        let mut bid = 0.;
        let iptt = i__ * 31 - (i__ + 1) * i__ / 2 + 1;
        for j in i__..=*ndeg / 2 {
            // OCCT: Getmmjcobi().plgcan[iptt + j + ideriv * 992 + 991]
            //   (0-based flat index over the 3968-entry plgcan table).
            bid += PLGCAN[(iptt + j + ideriv * 992 + 991) as usize] * poljac[j as usize];
        }
        polcan[(i__ * 2) as usize] = bid;
    }

    // --------------- Expression of terms of uneven degree ----------------
    if *ndeg == 0 {
        // OCCT goto L9999.
        if ibb >= 5 {
            sys_base::mgsomsg_("MMJACAN");
        }
        return;
    }

    for i__ in 0..=(*ndeg - 1) / 2 {
        let mut bid = 0.;
        let iptt = i__ * 31 - (i__ + 1) * i__ / 2 + 1;
        for j in i__..=(*ndeg - 1) / 2 {
            // OCCT: Getmmjcobi().plgcan[iptt + j + ((ideriv << 1) + 1) * 496
            //   + 991] (0-based flat index).
            bid += PLGCAN[(iptt + j + ((ideriv << 1) + 1) * 496 + 991) as usize]
                * poljac[(j + poljac_dim1) as usize];
        }
        polcan[((i__ << 1) + 1) as usize] = bid;
    }

    // -------------------------------- The end ----------------------------
    if ibb >= 5 {
        sys_base::mgsomsg_("MMJACAN");
    }
}

/// OCCT AdvApp2Var_MathBase::mmjaccv_ (AdvApp2Var_MathBase.cxx L5792-6136) -
/// passage from the normalized Jacobi base to the canonic base.
pub fn mmjaccv_(
    ncoef: &i32,
    ndim: &i32,
    ider: &i32,
    crvlgd: &[f64],
    polaux: &mut [f64],
    crvcan: &mut [f64],
) {
    // OCCT: polaux_dim1 = (*ncoef - 1) / 2 + 1 (polaux is a raw 0-based
    // workspace, no pointer adjustment).
    let polaux_dim1 = (*ncoef - 1) / 2 + 1;
    // OCCT: crvcan_dim1 = *ncoef; crvcan_offset = crvcan_dim1.
    let crvcan_dim1 = *ncoef;
    let crvcan_off = crvcan_dim1;
    // OCCT: crvlgd_dim1 = *ncoef; crvlgd_offset = crvlgd_dim1.
    let crvlgd_dim1 = *ncoef;
    let crvlgd_off = crvlgd_dim1;

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMJACCV ");
    }

    let mut ndeg = *ncoef - 1;

    for nd in 1..=*ndim {
        //  Loading of the auxiliary table.
        let mut ii = 0i32;
        for i__ in 0..=ndeg / 2 {
            // OCCT: polaux[i__] = crvlgd[ii + nd * crvlgd_dim1]
            //   (adjusted = raw[(nd - 1) * crvlgd_dim1 + ii]).
            polaux[i__ as usize] =
                crvlgd[((nd - 1) * crvlgd_dim1 + ii) as usize];
            ii += 2;
        }

        ii = 1;
        if ndeg >= 1 {
            for i__ in 0..=(ndeg - 1) / 2 {
                // OCCT: polaux[i__ + polaux_dim1] = crvlgd[ii + nd * dim1].
                polaux[(i__ + polaux_dim1) as usize] =
                    crvlgd[((nd - 1) * crvlgd_dim1 + ii) as usize];
                ii += 2;
            }
        }
        //   Call the routine of base change.
        // OCCT: mmjacan_(ider, &ndeg, polaux, &crvcan[nd * crvcan_dim1])
        //   (adjusted = raw[((nd - 1) * crvcan_dim1)..]).
        mmjacan_(ider, &mut ndeg, polaux, &mut crvcan[((nd - 1) * crvcan_dim1) as usize..]);
    }
}

/// OCCT AdvApp2Var_MathBase::mmmpocur_ (AdvApp2Var_MathBase.cxx L6596-6893) -
/// position of a point on curve (ncofmx, ndim).  `tabval` is 1-based
/// emulated.
pub fn mmmpocur_(
    ncofmx: &mut i32,
    ndim: &mut i32,
    ndeg: &mut i32,
    courbe: &[f64],
    tparam: &mut f64,
    tabval: &mut [f64],
) {
    // OCCT: courbe_dim1 = *ncofmx; courbe_offset = courbe_dim1 + 1.
    let courbe_dim1 = *ncofmx;
    let courbe_off = courbe_dim1 + 1;
    let c = |i: i32, nd: i32| -> usize {
        // OCCT: courbe[i + nd * courbe_dim1] - courbe_offset.
        ((i + nd * courbe_dim1) - courbe_off) as usize
    };

    if *ndeg < 1 {
        for nd in 1..=*ndim {
            // OCCT: tabval[nd] (1-based).
            tabval[(nd - 1) as usize] = 0.;
        }
    } else {
        for nd in 1..=*ndim {
            let mut fu = courbe[c(*ndeg, nd)];
            for i__ in (1..=*ndeg - 1).rev() {
                fu = fu * *tparam + courbe[c(i__, nd)];
            }
            tabval[(nd - 1) as usize] = fu;
        }
    }
}

/// OCCT mmtrpj0_ (AdvApp2Var_MathBase.cxx L8314-8422) - lowers the degree of
/// a Legendre-basis curve on (-1,1) with a given precision (stop at 1).
pub fn mmtrpj0_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    epsi3d: &mut f64,
    crvlgd: &[f64],
    ycvmax: &mut [f64],
    epstrc: &mut f64,
    ncfnew: &mut i32,
) {
    // OCCT: crvlgd_dim1 = *ncofmx; crvlgd_offset = crvlgd_dim1 + 1.
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    *ncfnew = 1;
    // ------------------- Init for error calculation ----------------------
    for i__ in 1..=*ndimen {
        ycvmax[(i__ - 1) as usize] = 0.;
    }
    *epstrc = 0.;
    let mut error;

    //   Cutting of coefficients.
    let ncut = 2;
    // ------ Loop on the series of Legendre :NCOEFF --> 2 (RBD) -----------
    for i__ in (ncut..=*ncoeff).rev() {
        //   Factor of renormalization.
        let mut bidon = ((i__ - 1) as f64 * 2. + 1.) / 2.;
        bidon = bidon.sqrt();
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvlgd[cl(i__, nd)].abs() * bidon;
        }
        //   Cutting is stopped if the norm becomes too great.
        error = mzsnorm_(ndimen, ycvmax);
        if error > *epsi3d {
            *ncfnew = i__;
            return; // OCCT goto L9999.
        }

        // ---  Max error cumulee when the I-th coeff is removed.
        *epstrc = error;
    }
}

/// OCCT mmtrpj2_ (AdvApp2Var_MathBase.cxx L8427-8595) - degree lowering in
/// the Jacobi basis of order 2 (ia = 2).  OCCT carries an inline `xmaxj[57]`
/// static identical to mma2jmx_'s XMAX2 (approxf2var_b.rs).
pub fn mmtrpj2_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    epsi3d: &mut f64,
    crvlgd: &[f64],
    ycvmax: &mut [f64],
    epstrc: &mut f64,
    ncfnew: &mut i32,
) {
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    //   Minimum degree that can be reached : Stop at IA (RBD). ------------
    let ia = 2;
    *ncfnew = ia;
    // Init for calculation of error.
    for i__ in 1..=*ndimen {
        ycvmax[(i__ - 1) as usize] = 0.;
    }
    *epstrc = 0.;
    let mut error;

    //   Cutting of coefficients.
    let ncut = ia + 1;
    // ------ Loop on the series of Jacobi :NCOEFF --> IA+1 (RBD) ----------
    'l300: for i__ in (ncut..=*ncoeff).rev() {
        //   Factor of renormalization.
        let bidon = XMAX2[(i__ - ncut) as usize];
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvlgd[cl(i__, nd)].abs() * bidon;
        }
        //   One stops to cut if the norm becomes too great.
        error = mzsnorm_(ndimen, ycvmax);
        if error > *epsi3d {
            *ncfnew = i__;
            break 'l300; // OCCT goto L400.
        }

        // --- Max error cumulated when the I-th coeff is removed.
        *epstrc = error;
    }

    // ------- Cutting of zero coeffs of interpolation (RBD) -------
    // OCCT L400:
    if *ncfnew == ia {
        let mut eps1 = 0.;
        mmeps1_(&mut eps1);
        for i__ in (2..=ia).rev() {
            let mut bid = 0.;
            for nd in 1..=*ndimen {
                bid += crvlgd[cl(i__, nd)].abs();
            }
            if bid > eps1 {
                *ncfnew = i__;
                return; // OCCT goto L9999.
            }
        }
        // --- If all coeffs can be removed, this is a point.
        *ncfnew = 1;
    }
}

/// OCCT mmtrpj4_ (AdvApp2Var_MathBase.cxx L8600-8767) - Jacobi basis of
/// order 4 (ia = 4).  OCCT carries an inline `xmaxj[55]` static identical to
/// mma2jmx_'s XMAX4 (approxf2var_b.rs).
pub fn mmtrpj4_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    epsi3d: &mut f64,
    crvlgd: &[f64],
    ycvmax: &mut [f64],
    epstrc: &mut f64,
    ncfnew: &mut i32,
) {
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    let ia = 4;
    *ncfnew = ia;
    for i__ in 1..=*ndimen {
        ycvmax[(i__ - 1) as usize] = 0.;
    }
    *epstrc = 0.;
    let mut error;

    let ncut = ia + 1;
    'l300: for i__ in (ncut..=*ncoeff).rev() {
        let bidon = XMAX4[(i__ - ncut) as usize];
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvlgd[cl(i__, nd)].abs() * bidon;
        }
        error = mzsnorm_(ndimen, ycvmax);
        if error > *epsi3d {
            *ncfnew = i__;
            break 'l300; // OCCT goto L400.
        }
        *epstrc = error;
    }

    // OCCT L400:
    if *ncfnew == ia {
        let mut eps1 = 0.;
        mmeps1_(&mut eps1);
        for i__ in (2..=ia).rev() {
            let mut bid = 0.;
            for nd in 1..=*ndimen {
                bid += crvlgd[cl(i__, nd)].abs();
            }
            if bid > eps1 {
                *ncfnew = i__;
                return; // OCCT goto L9999.
            }
        }
        *ncfnew = 1;
    }
}

/// OCCT mmtrpj6_ (AdvApp2Var_MathBase.cxx L8768-8941) - Jacobi basis of
/// order 6 (ia = 6).  OCCT carries an inline `xmaxj[53]` static identical to
/// mma2jmx_'s XMAX6 (approxf2var_b.rs).
pub fn mmtrpj6_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    epsi3d: &mut f64,
    crvlgd: &[f64],
    ycvmax: &mut [f64],
    epstrc: &mut f64,
    ncfnew: &mut i32,
) {
    let crvlgd_dim1 = *ncofmx;
    let crvlgd_off = crvlgd_dim1 + 1;
    let cl = |i: i32, nd: i32| -> usize {
        ((i + nd * crvlgd_dim1) - crvlgd_off) as usize
    };

    let ia = 6;
    *ncfnew = ia;
    for i__ in 1..=*ndimen {
        ycvmax[(i__ - 1) as usize] = 0.;
    }
    *epstrc = 0.;
    let mut error;

    let ncut = ia + 1;
    'l300: for i__ in (ncut..=*ncoeff).rev() {
        let bidon = XMAX6[(i__ - ncut) as usize];
        for nd in 1..=*ndimen {
            ycvmax[(nd - 1) as usize] += crvlgd[cl(i__, nd)].abs() * bidon;
        }
        error = mzsnorm_(ndimen, ycvmax);
        if error > *epsi3d {
            *ncfnew = i__;
            break 'l300; // OCCT goto L400.
        }
        *epstrc = error;
    }

    // OCCT L400:
    if *ncfnew == ia {
        let mut eps1 = 0.;
        mmeps1_(&mut eps1);
        for i__ in (2..=ia).rev() {
            let mut bid = 0.;
            for nd in 1..=*ndimen {
                bid += crvlgd[cl(i__, nd)].abs();
            }
            if bid > eps1 {
                *ncfnew = i__;
                return; // OCCT goto L9999.
            }
        }
        *ncfnew = 1;
    }
}

/// OCCT AdvApp2Var_MathBase::mmtrpjj_ (AdvApp2Var_MathBase.cxx L8943-9026) -
/// dispatch of the degree lowering on jord = 2*(IORDRE+1).
pub fn mmtrpjj_(
    ncofmx: &mut i32,
    ndimen: &mut i32,
    ncoeff: &mut i32,
    epsi3d: &mut f64,
    iordre: &mut i32,
    crvlgd: &[f64],
    ycvmax: &mut [f64],
    errmax: &mut f64,
    ncfnew: &mut i32,
) {
    let ia = (*iordre + 1) << 1;

    // OCCT: the crvlgd sub-pointer &crvlgd[crvlgd_offset] resolves to the
    // raw address (the rcad slice is the physical buffer); ycvmax is passed
    // as &ycvmax[1] (1-based emulated, raw index 0).
    if ia == 0 {
        mmtrpj0_(ncofmx, ndimen, ncoeff, epsi3d, crvlgd, ycvmax, errmax, ncfnew);
    } else if ia == 2 {
        mmtrpj2_(ncofmx, ndimen, ncoeff, epsi3d, crvlgd, ycvmax, errmax, ncfnew);
    } else if ia == 4 {
        mmtrpj4_(ncofmx, ndimen, ncoeff, epsi3d, crvlgd, ycvmax, errmax, ncfnew);
    } else {
        mmtrpj6_(ncofmx, ndimen, ncoeff, epsi3d, crvlgd, ycvmax, errmax, ncfnew);
    }
}
