//! OCCT AdvApp2Var_ApproxF2var (AdvApp2Var_ApproxF2var.cxx) - part F: the
//! function-discretization engines (mma2ds1_, mma2ds2_) and the
//! cut-managed approximation driver mma2fnc_ consumed by AdvApp2Var_Iso.
//!
//! Encoding note: as in parts A-E, the OCCT Fortran 1-based index formulas
//! are kept on 0-based slices by subtracting the OCCT virtual origin inside
//! the flat-index expression; workspace zones are expressed with
//! split_at_mut over the raw allocation.

use super::approxf2var_c::mma1cdi_;
use super::approxf2var_c::mma1cnt_;
use super::approxf2var_c::mma1fdi_;
use super::approxf2var_c::mma1fer_;
use super::approxf2var_c::mma1jak_;
use super::approxf2var_c::mma1noc_;
use super::approxf2var_c::mma1nop_;
use super::approxf2var_c::EvaluatorFunc2Var;
use super::math_base_b::mmfmtb1_;
use super::math_base_b::mmjacan_;
use super::math_base_b::mmveps3_;
use super::math_base_b::mmfmca8_;
use super::math_base_b::mmapcmp_;
use super::sys_base;

// ---------------------------------------------------------------------------
// mma2ds2_ (AdvApp2Var_ApproxF2var.cxx L5789-6207)
// ---------------------------------------------------------------------------

/// OCCT mma2ds2_ (AdvApp2Var_ApproxF2var.cxx L5789-6207) - discretization of
/// F(u,v) on the Legendre roots of degree NBPNTU by U with Vj fixed.
/// `uintfn` / `vintfn` / `ttable` / `urootb` / `vrootb` are 1-based
/// emulated; the sosotb-family buffers and `fpntab` are raw physical
/// slices.
#[allow(clippy::too_many_arguments)]
pub fn mma2ds2_(
    ndimen: &mut i32,
    uintfn: &[f64],
    vintfn: &[f64],
    foncnp: &dyn EvaluatorFunc2Var,
    nbpntu: &mut i32,
    nbpntv: &mut i32,
    urootb: &[f64],
    vrootb: &[f64],
    iiuouv: &mut i32,
    sosotb: &mut [f64],
    disotb: &mut [f64],
    soditb: &mut [f64],
    diditb: &mut [f64],
    fpntab: &mut [f64],
    ttable: &mut [f64],
    iercod: &mut i32,
) {
    let c__0: i32 = 0;

    // OCCT: sosotb_dim1 = *nbpntu / 2 + 1; sosotb_dim2 = *nbpntv / 2 + 1.
    let sosotb_dim1 = (*nbpntu / 2 + 1) as i32;
    let sosotb_dim2 = (*nbpntv / 2 + 1) as i32;
    let diditb_dim1 = (*nbpntu / 2 + 1) as i32;
    let diditb_dim2 = (*nbpntv / 2 + 1) as i32;
    let soditb_dim1 = (*nbpntu / 2) as i32;
    let soditb_dim2 = (*nbpntv / 2) as i32;
    // OCCT: fpntab_dim1 = *ndimen (fpntab_offset = dim1 + 1).
    let fpntab_dim1 = *ndimen;

    let fp = |id: i32, col: i32| -> usize {
        // OCCT: fpntab[id + col * fpntab_dim1] - (dim1 + 1).
        ((id - 1) + (col - 1) * fpntab_dim1) as usize
    };
    let ss = |iu: i32, jj: i32, id: i32| -> usize {
        // OCCT: sosotb[iu + (jj + id * dim2) * dim1] - dim1 * dim2.
        ((iu - 1) + ((jj - 1) + (id - 1) * sosotb_dim2) * sosotb_dim1) as usize
    };
    let dd = |iu: i32, jj: i32, id: i32| -> usize {
        ((iu - 1) + ((jj - 1) + (id - 1) * diditb_dim2) * diditb_dim1) as usize
    };
    let sd = |iu: i32, jj: i32, id: i32| -> usize {
        // OCCT: soditb/disotb[iu + (jj + id * dim2) * dim1]
        //   - dim1 * (dim2 + 1) - 1.
        ((iu - 1) + ((jj - 1) + (id - 1) * soditb_dim2) * soditb_dim1) as usize
    };
    // Slot-0 (iu = 0) forms.
    let ss0 = |jj: i32, id: i32| -> usize {
        // OCCT: sosotb[(jj + id * dim2) * dim1] - dim1 * dim2.
        ((jj - 1) + (id - 1) * sosotb_dim2) as usize * sosotb_dim1 as usize
    };
    let dd0 = |jj: i32, id: i32| -> usize {
        ((jj - 1) + (id - 1) * diditb_dim2) as usize * diditb_dim1 as usize
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2DS2");
    }
    *iercod = 0;

    // OCCT: uintfn/vintfn are 1-based emulated (raw [0]/[1]).
    let alinu = (uintfn[1] - uintfn[0]) / 2.;
    let blinu = (uintfn[1] + uintfn[0]) / 2.;
    let alinv = (vintfn[1] - vintfn[0]) / 2.;
    let blinv = (vintfn[1] + vintfn[0]) / 2.;

    let mut dbfn1 = [0.0f64; 2];
    let mut dbfn2 = [0.0f64; 2];
    if *iiuouv == 1 {
        dbfn1[0] = vintfn[0];
        dbfn1[1] = vintfn[1];
        dbfn2[0] = uintfn[0];
        dbfn2[1] = uintfn[1];
    } else {
        dbfn1[0] = uintfn[0];
        dbfn1[1] = uintfn[1];
        dbfn2[0] = vintfn[0];
        dbfn2[1] = vintfn[1];
    }

    // -------- Discretization by U on the roots of Legendre polynom --------
    // ---------------- of degree NBPNTU, with Vj fixed  --------------------
    let nuroo = *nbpntu / 2;
    let nvroo = *nbpntv / 2;
    let jdec = (*nbpntu + 1) / 2;

    // ----------- Loading of parameters of discretization by U -------------
    for iu in 1..=*nbpntu {
        // OCCT: ttable[iu] (1-based emulated).
        ttable[(iu - 1) as usize] = blinu + alinu * urootb[(iu - 1) as usize];
    }

    // -------------- For Vj fixed, negative root of Legendre -------------
    for iv in 1..=nvroo {
        let mut tcons = blinv + alinv * vrootb[(iv - 1) as usize];
        // OCCT: foncnp.Evaluate(ndimen, dbfn1, dbfn2, iiuouv, &tcons,
        //   nbpntu, &ttable[1], &c__0, &c__0,
        //   &fpntab[fpntab_offset], iercod) (sub-pointer raw 0).
        foncnp.evaluate(
            ndimen,
            &dbfn1,
            &dbfn2,
            iiuouv,
            &tcons,
            nbpntu,
            ttable,
            &c__0,
            &c__0,
            fpntab,
            iercod,
        );
        if *iercod > 0 {
            finish_mma2ds2_(ldbg, iercod);
            return;
        }
        for id in 1..=*ndimen {
            for iu in 1..=nuroo {
                let up = fpntab[fp(id, iu + jdec)];
                let um = fpntab[fp(id, nuroo - iu + 1)];
                let jj = nvroo - iv + 1;
                sosotb[ss(iu, jj, id)] += up + um;
                disotb[sd(iu, jj, id)] += up - um;
                soditb[sd(iu, jj, id)] -= up + um;
                diditb[dd(iu, jj, id)] -= up - um;
            }
            if *nbpntu % 2 != 0 {
                let up = fpntab[fp(id, jdec)];
                let jj = nvroo - iv + 1;
                sosotb[ss0(jj, id)] += up;
                diditb[dd0(jj, id)] -= up;
            }
        }
    }

    // --------- For Vj = 0 (uneven NBPNTV), discretization by U -----------
    if *nbpntv % 2 != 0 {
        let mut tcons = blinv;
        foncnp.evaluate(
            ndimen,
            &dbfn1,
            &dbfn2,
            iiuouv,
            &tcons,
            nbpntu,
            ttable,
            &c__0,
            &c__0,
            fpntab,
            iercod,
        );
        if *iercod > 0 {
            finish_mma2ds2_(ldbg, iercod);
            return;
        }
        for id in 1..=*ndimen {
            for iu in 1..=nuroo {
                let up = fpntab[fp(id, jdec + iu)];
                let um = fpntab[fp(id, nuroo - iu + 1)];
                sosotb[ss(iu, 0, id)] += up + um;
                diditb[dd(iu, 0, id)] += up - um;
            }
            if *nbpntu % 2 != 0 {
                let up = fpntab[fp(id, jdec)];
                sosotb[ss0(0, id)] += up;
            }
        }
    }

    // -------------- For Vj fixed, positive root of Legendre -------------
    for iv in 1..=nvroo {
        let mut tcons = alinv * vrootb[(((*nbpntv + 1) / 2 + iv) - 1) as usize] + blinv;
        foncnp.evaluate(
            ndimen,
            &dbfn1,
            &dbfn2,
            iiuouv,
            &tcons,
            nbpntu,
            ttable,
            &c__0,
            &c__0,
            fpntab,
            iercod,
        );
        if *iercod > 0 {
            finish_mma2ds2_(ldbg, iercod);
            return;
        }
        for id in 1..=*ndimen {
            for iu in 1..=nuroo {
                let up = fpntab[fp(id, iu + jdec)];
                let um = fpntab[fp(id, nuroo - iu + 1)];
                sosotb[ss(iu, iv, id)] += up + um;
                disotb[sd(iu, iv, id)] += up - um;
                soditb[sd(iu, iv, id)] += up + um;
                diditb[dd(iu, iv, id)] += up - um;
            }
            if *nbpntu % 2 != 0 {
                let up = fpntab[fp(id, jdec)];
                sosotb[ss0(iv, id)] += up;
                diditb[dd0(iv, id)] += up;
            }
        }
    }

    // ------------------------------ The end -------------------------------
    finish_mma2ds2_(ldbg, iercod);
}

/// OCCT mma2ds2_ L9999 tail (iercod += 100 + maermsg_ + mgsomsg_).
fn finish_mma2ds2_(ldbg: bool, iercod: &mut i32) {
    if *iercod > 0 {
        *iercod += 100;
        sys_base::maermsg_("MMA2DS2", iercod);
    }
    if ldbg {
        sys_base::mgsomsg_("MMA2DS2");
    }
}

// ---------------------------------------------------------------------------
// mma2ds1_ (AdvApp2Var_ApproxF2var.cxx L5383-5787)
// ---------------------------------------------------------------------------

/// OCCT mma2ds1_ (AdvApp2Var_ApproxF2var.cxx L5383-5787) - discretization of
/// F(u,v) on the roots of Legendre polynoms, dispatching on the iso type.
#[allow(clippy::too_many_arguments)]
pub fn mma2ds1_(
    ndimen: &mut i32,
    uintfn: &[f64],
    vintfn: &[f64],
    foncnp: &dyn EvaluatorFunc2Var,
    nbpntu: &mut i32,
    nbpntv: &mut i32,
    urootb: &[f64],
    vrootb: &[f64],
    isofav: &mut i32,
    sosotb: &mut [f64],
    disotb: &mut [f64],
    soditb: &mut [f64],
    diditb: &mut [f64],
    fpntab: &mut [f64],
    ttable: &mut [f64],
    iercod: &mut i32,
) {
    let sosotb_dim1 = (*nbpntu / 2 + 1) as i32;
    let sosotb_dim2 = (*nbpntv / 2 + 1) as i32;
    let diditb_dim1 = (*nbpntu / 2 + 1) as i32;
    let diditb_dim2 = (*nbpntv / 2 + 1) as i32;
    let soditb_dim1 = (*nbpntu / 2) as i32;
    let soditb_dim2 = (*nbpntv / 2) as i32;

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2DS1");
    }
    *iercod = 0;
    let iuouv = if *isofav < 1 || *isofav > 2 {
        2
    } else {
        *isofav
    };
    let mut iuouv = iuouv;

    // --------- Discretization by U on the roots of the polynom of --------
    // --------------- Legendre of degree NBPNTU, iso-V by iso-V -----------
    if iuouv == 2 {
        mma2ds2_(
            ndimen,
            uintfn,
            vintfn,
            foncnp,
            nbpntu,
            nbpntv,
            urootb,
            vrootb,
            &mut iuouv,
            sosotb,
            disotb,
            soditb,
            diditb,
            fpntab,
            ttable,
            iercod,
        );
    } else {
        // --> Inversion of indices of tables (in-place transposes).
        let mut ibid1 = 0;
        let mut ibid2 = 0;
        for nd in 1..=*ndimen {
            let isz1 = *nbpntu / 2 + 1;
            let isz2 = *nbpntv / 2 + 1;
            let base = ((nd - 1) * sosotb_dim2 * sosotb_dim1) as usize;
            // OCCT: mmfmtb1_(&isz1, &sosotb[nd*dim2*dim1], &isz1, &isz2,
            //   &isz2, &sosotb[nd*dim2*dim1], &ibid1, &ibid2, iercod)
            //   (in-place: table1 == table2; the callee stages through its
            //   WORK block, the rcad encoding snapshots the input;
            //   maxsz1 aliases isize1 and maxsz2 aliases jsize1, so the
            //   rcad encoding splits them into locals).
            let snap: Vec<f64> = sosotb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut sosotb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
            let base = ((nd - 1) * diditb_dim2 * diditb_dim1) as usize;
            let snap: Vec<f64> = diditb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut diditb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
            let isz1 = *nbpntu / 2;
            let isz2 = *nbpntv / 2;
            // OCCT: &soditb[(nd * dim2 + 1) * dim1 + 1] resolves to the raw
            // address (nd - 1) * dim1 * dim2.
            let base = ((nd - 1) * soditb_dim2 * soditb_dim1) as usize;
            let snap: Vec<f64> = soditb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut soditb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
            let base = ((nd - 1) * soditb_dim2 * soditb_dim1) as usize;
            let snap: Vec<f64> = disotb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut disotb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
        }

        mma2ds2_(
            ndimen,
            vintfn,
            uintfn,
            foncnp,
            nbpntv,
            nbpntu,
            vrootb,
            urootb,
            &mut iuouv,
            sosotb,
            soditb,
            disotb,
            diditb,
            fpntab,
            ttable,
            iercod,
        );

        // --> Inversion of indices of tables
        for nd in 1..=*ndimen {
            let isz1 = *nbpntv / 2 + 1;
            let isz2 = *nbpntu / 2 + 1;
            let base = ((nd - 1) * sosotb_dim2 * sosotb_dim1) as usize;
            let snap: Vec<f64> = sosotb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut sosotb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
            let base = ((nd - 1) * diditb_dim2 * diditb_dim1) as usize;
            let snap: Vec<f64> = diditb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut diditb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
            let isz1 = *nbpntv / 2;
            let isz2 = *nbpntu / 2;
            let base = ((nd - 1) * soditb_dim2 * soditb_dim1) as usize;
            let snap: Vec<f64> = soditb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut soditb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
            let base = ((nd - 1) * soditb_dim2 * soditb_dim1) as usize;
            let snap: Vec<f64> = disotb[base..base + (isz1 * isz2) as usize].to_vec();
            ds1_transpose_(isz1, isz2, &snap, &mut disotb[base..], &mut ibid1, &mut ibid2, iercod);
            if *iercod > 0 {
                finish_mma2ds1_(ldbg, iercod);
                return;
            }
        }
    }

    // ------------------------------ The end -------------------------------
    finish_mma2ds1_(ldbg, iercod);
}

/// OCCT mma2ds1_ L9999 tail.
fn finish_mma2ds1_(ldbg: bool, iercod: &mut i32) {
    if *iercod > 0 {
        *iercod += 100;
        sys_base::maermsg_("MMA2DS1", iercod);
    }
    if ldbg {
        sys_base::mgsomsg_("MMA2DS1");
    }
}

/// Encoding helper for the OCCT in-place mmfmtb1_ calls of mma2ds1_
/// (maxsz1 aliases isize1 and maxsz2 aliases jsize1 - the callee only
/// reads those, so the rcad encoding splits them into locals).
#[allow(clippy::too_many_arguments)]
fn ds1_transpose_(
    isz1: i32,
    isz2: i32,
    snap: &[f64],
    target: &mut [f64],
    ibid1: &mut i32,
    ibid2: &mut i32,
    iercod: &mut i32,
) {
    let mut a1 = isz1;
    let mut a2 = isz1;
    let mut b1 = isz2;
    let mut b2 = isz2;
    mmfmtb1_(&mut a1, snap, &mut a2, &mut b1, &mut b2, target, ibid1, ibid2, iercod);
}

// ---------------------------------------------------------------------------
// mma2fnc_ (AdvApp2Var_ApproxF2var.cxx L6624-7283)
// ---------------------------------------------------------------------------

/// OCCT mma2fnc_ (AdvApp2Var_ApproxF2var.cxx L6624-7283) - approximation
/// with cuts of F(u,v) by polynomial curves.  `uvfonc` / `rootlg` / `epsapr`
/// / `ndimse` / `ncoeff` are 1-based emulated; somtab / diftab / contr1 /
/// contr2 / tabdec / errmax / errmoy / courbe are raw physical slices.
#[allow(clippy::too_many_arguments)]
pub fn mma2fnc_(
    ndimen: &mut i32,
    nbsesp: &mut i32,
    ndimse: &[i32],
    uvfonc: &[f64],
    foncnp: &dyn EvaluatorFunc2Var,
    tconst: &mut f64,
    isofav: &mut i32,
    nbroot: &mut i32,
    rootlg: &[f64],
    iordre: &mut i32,
    ideriv: &mut i32,
    ndgjac: &mut i32,
    nbcrmx: &mut i32,
    ncflim: &mut i32,
    epsapr: &[f64],
    ncoeff: &mut [i32],
    courbe: &mut [f64],
    nbcrbe: &mut i32,
    somtab: &mut [f64],
    diftab: &mut [f64],
    contr1: &mut [f64],
    contr2: &mut [f64],
    tabdec: &mut [f64],
    errmax: &mut [f64],
    errmoy: &mut [f64],
    iercod: &mut i32,
) {
    let c__8: i32 = 8;

    // OCCT: courbe_dim1 = *ncflim; courbe_dim2 = *ndimen;
    //       courbe_offset = dim1 * (dim2 + 1) + 1.
    let courbe_dim1 = *ncflim;
    let courbe_dim2 = *ndimen;
    // OCCT: somtab_dim1 = *nbroot / 2 + 1; somtab_dim2 = *ndimen.
    let somtab_dim1 = (*nbroot / 2 + 1) as i32;
    let somtab_dim2 = *ndimen;
    let diftab_dim1 = somtab_dim1;
    let diftab_dim2 = *ndimen;
    // OCCT: contrX_dim1 = *ndimen; contrX_dim2 = *iordre + 2.
    let contr1_dim1 = *ndimen;
    let contr1_dim2 = *iordre + 2;
    let contr2_dim1 = *ndimen;
    let contr2_dim2 = *iordre + 2;
    // OCCT: errmax_dim1 = *nbsesp.
    let errmax_dim1 = *nbsesp;

    let ibb = sys_base::mnfndeb_();
    if ibb >= 1 {
        sys_base::mgenmsg_("MMA2FNC");
    }
    *iercod = 0;

    // ---------------- Set to zero the coefficients of CURBE --------------
    let ilong = *ndimen * *ncflim * *nbcrmx;
    sys_base::mvriraz_(ilong, courbe);

    // -------------------------- Checking of entries ----------------------
    let mut eps3 = 0.;
    mmveps3_(&mut eps3);
    // OCCT: uvfonc[4] - uvfonc[3] / uvfonc[6] - uvfonc[5] (raw [1]/[3]).
    if (uvfonc[1] - uvfonc[0]).abs() < eps3 {
        // OCCT goto L9100.
        *iercod = 1;
        finish_mma2fnc_(ibb, iercod);
        return;
    }
    if (uvfonc[3] - uvfonc[2]).abs() < eps3 {
        *iercod = 1;
        finish_mma2fnc_(ibb, iercod);
        return;
    }

    let mut uv11 = [-1., 1., -1., 1.0f64];

    // ------------- Preparation of parameters of discretisation -----------
    // -- Allocation of a table of parameters and points of discretisation --
    let isz1 = *nbroot + 2;
    let mut ibid1 = *ndimen * (*nbroot + 2);
    let ibid2 = ((*iordre + 1) << 1) * *nbroot;
    let mut isz2 = ibid1.max(ibid2);
    ibid1 = (((*ncflim - 1) / 2 + 1) << 1) * *ndimen;
    isz2 = ibid1.max(isz2);
    // --> To return the polynoms of hermit.
    let isz3 = ((*iordre + 1) << 2) * (*iordre + 1);
    // --> For the Gauss  coeff. of integration.
    let isz4 = (*nbroot / 2 + 1) * (*ndgjac + 1 - ((*iordre + 1) << 1));
    // --> For the coeff of the curve in the base of Jacobi
    let isz5 = (*ndgjac + 1) * *ndimen;

    let mut ndwrk = isz1 + isz2 + isz3 + isz4 + isz5;
    let mut wrkar: Vec<f64> = Vec::new();
    let mut iofwr: isize = 0;
    let mut ier = 0;
    let mut iunit = c__8;
    sys_base::mcrrqst_(&mut iunit, &mut ndwrk, &mut wrkar, &mut iofwr, &mut ier);
    if ier > 0 {
        // OCCT goto L9013.
        *iercod = 13;
        finish_mma2fnc_(ibb, iercod);
        return;
    }
    let ipt1 = isz1;
    let ipt2 = ipt1 + isz2;
    let ipt3 = ipt2 + isz3;
    let ipt4 = ipt3 + isz4;

    // ------------------ Initialisation of management of cuts ---------
    let mut uvpav = [0.0f64; 4];
    if *isofav == 1 {
        uvpav[0] = uvfonc[0];
        uvpav[1] = uvfonc[1];
        tabdec[0] = uvfonc[2];
        tabdec[1] = uvfonc[3];
    } else if *isofav == 2 {
        tabdec[0] = uvfonc[0];
        tabdec[1] = uvfonc[1];
        uvpav[2] = uvfonc[2];
        uvpav[3] = uvfonc[3];
    } else {
        *iercod = 1;
        finish_mma2fnc_(ibb, iercod);
        return;
    }

    let mut nupil: i32 = 1;
    *nbcrbe = 0;

    //                       APPROXIMATION WITH CUTS
    // OCCT L1000 loop:
    loop {
        // --> When the top is reached, this is the end !
        if nupil - *nbcrbe == 0 {
            // OCCT goto L9900.
            break;
        }
        let ncb1 = *nbcrbe + 1;
        if *isofav == 1 {
            uvpav[2] = tabdec[*nbcrbe as usize];
            uvpav[3] = tabdec[(*nbcrbe + 1) as usize];
        } else if *isofav == 2 {
            uvpav[0] = tabdec[*nbcrbe as usize];
            uvpav[1] = tabdec[(*nbcrbe + 1) as usize];
        } else {
            *iercod = 1;
            break;
        }

        // -------------------- Normalization of parameters --------------------
        // OCCT: mma1nop_(nbroot, &rootlg[1], uvpav, isofav, wrkar_off, &ier)
        //   (ttable = the wrkar base zone).
        {
            let (w_ttable, _) = wrkar.split_at_mut(isz1 as usize);
            mma1nop_(nbroot, rootlg, &uvpav, isofav, w_ttable, &mut ier);
        }
        if ier > 0 {
            // OCCT goto L9100.
            *iercod = 1;
            break;
        }

        // -------------------- Discretisation of FONCNP ------------------------
        // Sub-buffer raw bases for ncb1.
        let somtab_base = ((ncb1 - 1) * somtab_dim1 * somtab_dim2) as usize;
        let diftab_base = ((ncb1 - 1) * diftab_dim1 * diftab_dim2) as usize;
        let contr1_base = ((ncb1 - 1) * contr1_dim1 * contr1_dim2) as usize;
        let contr2_base = ((ncb1 - 1) * contr2_dim1 * contr2_dim2) as usize;
        {
            let (w_ttable, w_rest) = wrkar.split_at_mut(ipt1 as usize);
            let (w_fpntab, _) = w_rest.split_at_mut(ipt2 as usize - ipt1 as usize);
            mma1fdi_(
                ndimen,
                &uvpav,
                foncnp,
                isofav,
                tconst,
                nbroot,
                w_ttable,
                iordre,
                ideriv,
                w_fpntab,
                &mut somtab[somtab_base..],
                &mut diftab[diftab_base..],
                &mut contr1[contr1_base..],
                &mut contr2[contr2_base..],
                iercod,
            );
        }
        if *iercod > 0 {
            // OCCT goto L9900.
            break;
        }

        // ----------- Cut the discretisation of constraints ------------
        if *iordre >= 0 {
            // OCCT: mma1cdi_(ndimen, nbroot, &rootlg[1], iordre,
            //   &contr1[(ncb1*dim2+1)*dim1+1], &contr2[...],
            //   &somtab[...], &diftab[...], &wrkar_off[ipt1],
            //   &wrkar_off[ipt2], &ier) (sub-pointers raw addresses).
            let (w_pre_ipt2, w_rest) = wrkar.split_at_mut(ipt2 as usize);
            let (_, w_fpntab) = w_pre_ipt2.split_at_mut(ipt1 as usize);
            mma1cdi_(
                ndimen,
                nbroot,
                rootlg,
                iordre,
                &mut contr1[contr1_base..],
                &mut contr2[contr2_base..],
                &mut somtab[somtab_base..],
                &mut diftab[diftab_base..],
                w_fpntab,
                w_rest,
                &mut ier,
            );
            if ier > 0 {
                // OCCT goto L9100.
                *iercod = 1;
                break;
            }
        }

        // -------------------- Calculate the curve of approximation -----------
        // OCCT: mma1jak_(ndimen, nbroot, iordre, ndgjac,
        //   &somtab[(ncb1*dim2+1)*dim1], &diftab[...], &wrkar_off[ipt3],
        //   &wrkar_off[ipt4], &ier).
        {
            let (_, w_rest) = wrkar.split_at_mut(ipt3 as usize);
            let (w_cgauss, w_crvjac) = w_rest.split_at_mut(ipt4 as usize - ipt3 as usize);
            let somtab_sub = &somtab[somtab_base..];
            let diftab_sub = &diftab[diftab_base..];
            mma1jak_(
                ndimen,
                nbroot,
                iordre,
                ndgjac,
                somtab_sub,
                diftab_sub,
                w_cgauss,
                w_crvjac,
                &mut ier,
            );
        }
        if ier > 0 {
            // OCCT goto L9100.
            *iercod = 1;
            break;
        }

        // ---------------- Add polynom of interpolation -------------------
        if *iordre >= 0 {
            // OCCT: mma1cnt_(ndimen, iordre, contr1_sub, contr2_sub,
            //   &wrkar_off[ipt2], ndgjac, &wrkar_off[ipt4]).
            let (a, b) = wrkar.split_at_mut(ipt3 as usize);
            let hermit = &mut a[ipt2 as usize..];
            let crvjac = &mut b[..isz5 as usize];
            mma1cnt_(
                ndimen,
                iordre,
                &contr1[contr1_base..],
                &contr2[contr2_base..],
                hermit,
                ndgjac,
                crvjac,
            );
        }

        // --------------- Calculate Max and Average error ---------------------
        // OCCT: mma1fer_(ndimen, nbsesp, &ndimse[1], iordre, ndgjac,
        //   &wrkar_off[ipt4], ncflim, &epsapr[1], &wrkar_off[ipt1],
        //   &errmax[ncb1*dim1+1], &errmoy[ncb1*dim1+1], &ncoeff[ncb1], &ier).
        {
            let (a, b) = wrkar.split_at_mut(ipt4 as usize);
            let ycvmax = &mut a[ipt1 as usize..];
            let crvjac = &mut b[..isz5 as usize];
            let emax_base = ((ncb1 - 1) * errmax_dim1) as usize;
            let mut ncoeff_ncb1 = ncoeff[(ncb1 - 1) as usize];
            mma1fer_(
                ndimen,
                nbsesp,
                ndimse,
                iordre,
                ndgjac,
                crvjac,
                ncflim,
                epsapr,
                ycvmax,
                &mut errmax[emax_base..],
                &mut errmoy[emax_base..],
                &mut ncoeff_ncb1,
                &mut ier,
            );
            ncoeff[(ncb1 - 1) as usize] = ncoeff_ncb1;
        }
        if ier > 0 {
            // OCCT goto L9100.
            *iercod = 1;
            break;
        }

        if ier == 0 || (ier == -1 && nupil == *nbcrmx) {
            // ----------------------- Result compression -------------------
            if ier == -1 {
                *iercod = -1;
            }
            let mut ncfja = *ndgjac + 1;
            // -> Compression of result in WRKAR(IPT1)
            // OCCT: mmapcmp_(ndimen, &ncfja, &ncoeff[ncb1],
            //   &wrkar_off[ipt4], &wrkar_off[ipt1]).
            {
                let (a, b) = wrkar.split_at_mut(ipt2 as usize);
                let crvnew = &mut a[ipt1 as usize..];
                let crvold = &mut b[..isz5 as usize];
                let mut ncoeff_ncb1 = ncoeff[(ncb1 - 1) as usize];
                mmapcmp_(ndimen, &mut ncfja, &mut ncoeff_ncb1, crvold, crvnew);
                ncoeff[(ncb1 - 1) as usize] = ncoeff_ncb1;
            }
            let ilong = *ndimen * *ncflim;
            {
                let (_, b) = wrkar.split_at_mut(ipt2 as usize);
                let crvold = &mut b[..ilong as usize];
                sys_base::mvriraz_(ilong, crvold);
            }
            // -> Passage to canonic base (-1,1) (result in WRKAR(IPT4)).
            let mut ndgre = ncoeff[(ncb1 - 1) as usize] - 1;
            for nd in 1..=*ndimen {
                let iptt = ipt1 + ((nd - 1) << 1) * (ndgre / 2 + 1);
                let jptt = ipt4 + (nd - 1) * ncoeff[(ncb1 - 1) as usize];
                // OCCT: mmjacan_(iordre, &ndgre, &wrkar_off[iptt],
                //   &wrkar_off[jptt]).
                let (a, b) = wrkar.split_at_mut(ipt4 as usize);
                let poljac = &mut a[iptt as usize..];
                let polcan = &mut b[jptt as usize - ipt4 as usize..];
                mmjacan_(iordre, &mut ndgre, poljac, polcan);
            }

            // -> Store the calculated curve
            // OCCT: mmfmca8_(&ncoeff[ncb1], ndimen, &ibid1, ncflim, ndimen,
            //   &ibid1, &wrkar_off[ipt4], &courbe[ncb1 sub]) with ibid1 = 1.
            let ibid1 = 1;
            {
                let (_, b) = wrkar.split_at_mut(ipt2 as usize);
                let tabini = &b[..(ncoeff[(ncb1 - 1) as usize] * *ndimen) as usize];
                let courbe_base = ((ncb1 - 1) * courbe_dim1 * courbe_dim2) as usize;
                let ncoeff_ncb1 = ncoeff[(ncb1 - 1) as usize];
                mmfmca8_(
                    &ncoeff_ncb1,
                    ndimen,
                    &ibid1,
                    ncflim,
                    ndimen,
                    &ibid1,
                    tabini,
                    &mut courbe[courbe_base..],
                );
            }

            // -> Before normalization of constraints on (-1,1), recalculate
            //    the true constraints.
            for ii in 0..=*iordre {
                // OCCT: mma1noc_(uv11, ndimen, &ii, contr1_sub, uvpav,
                //   isofav, ideriv, contr1_sub) (in/out same buffer).
                let base = ((ii + (ncb1 - 1) * contr1_dim2) * contr1_dim1) as usize;
                let snap: Vec<f64> = contr1[base..].to_vec();
                let mut ii_ = ii;
                mma1noc_(
                    &uv11,
                    ndimen,
                    &mut ii_,
                    &snap,
                    &uvpav,
                    isofav,
                    ideriv,
                    &mut contr1[base..],
                );
                let base = ((ii + (ncb1 - 1) * contr2_dim2) * contr2_dim1) as usize;
                let snap: Vec<f64> = contr2[base..].to_vec();
                let mut ii_ = ii;
                mma1noc_(
                    &uv11,
                    ndimen,
                    &mut ii_,
                    &snap,
                    &uvpav,
                    isofav,
                    ideriv,
                    &mut contr2[base..],
                );
            }
            let mut ii = 0;
            let ibid1 = (*nbroot / 2 + 1) * *ndimen;
            {
                let base = somtab_base;
                let snap: Vec<f64> = somtab[base..base + ibid1 as usize].to_vec();
                mma1noc_(
                    &uv11,
                    &ibid1,
                    &ii,
                    &snap,
                    &uvpav,
                    isofav,
                    ideriv,
                    &mut somtab[base..],
                );
                let base = diftab_base;
                let snap: Vec<f64> = diftab[base..base + ibid1 as usize].to_vec();
                mma1noc_(
                    &uv11,
                    &ibid1,
                    &ii,
                    &snap,
                    &uvpav,
                    isofav,
                    ideriv,
                    &mut diftab[base..],
                );
            }
            ii = 0;
            for nd in 1..=*ndimen {
                let base =
                    (((nd - 1) + (ncb1 - 1) * courbe_dim2) * courbe_dim1) as usize;
                let snap: Vec<f64> = courbe[base..].to_vec();
                let ncoeff_ncb1 = ncoeff[(ncb1 - 1) as usize];
                mma1noc_(
                    &uv11,
                    &ncoeff_ncb1,
                    &ii,
                    &snap,
                    &uvpav,
                    isofav,
                    ideriv,
                    &mut courbe[base..],
                );
            }

            // -> Update the nb of already created curves
            *nbcrbe += 1;

            // -> ...otherwise try to cut the current interval in 2...
        } else {
            let tmil = (tabdec[*nbcrbe as usize + 1] + tabdec[*nbcrbe as usize]) / 2.;
            let ideb = *nbcrbe + 1;
            let ideb1 = ideb + 1;
            let ilong = (nupil - *nbcrbe) << 3;
            // OCCT: mcrfill_(&ilong, &tabdec[ideb], &tabdec[ideb1]) -
            // overlapping forward shift; the rcad encoding stages the
            // source (mcrfill_ is memmove-safe in OCCT).
            let n = (ilong / 8) as usize;
            let src = (ideb - 1) as usize;
            let tmp: Vec<f64> = tabdec[src..src + n].to_vec();
            tabdec[ideb1 as usize - 1..ideb1 as usize - 1 + n].copy_from_slice(&tmp);
            tabdec[(ideb - 1) as usize] = tmil;
            nupil += 1;
        }
    }

    // OCCT L9900:
    {
        let mut ier2 = 0;
        let mut isize_ = 0;
        let mut iofst: isize = 0;
        let mut t: Vec<f64> = Vec::new();
        let mut iunit = c__8;
        sys_base::mcrdelt_(&mut iunit, &mut isize_, &mut t, &mut iofst, &mut ier2);
        if ier2 > 0 {
            *iercod = 13;
        }
    }
    finish_mma2fnc_(ibb, iercod);
}

/// OCCT mma2fnc_ L9999 tail.
fn finish_mma2fnc_(ibb: i32, iercod: &mut i32) {
    if *iercod != 0 {
        sys_base::maermsg_("MMA2FNC", iercod);
    }
    if ibb >= 2 {
        sys_base::mgsomsg_("MMA2FNC");
    }
}
