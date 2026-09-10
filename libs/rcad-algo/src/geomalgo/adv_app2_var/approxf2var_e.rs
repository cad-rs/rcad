//! OCCT AdvApp2Var_ApproxF2var (AdvApp2Var_ApproxF2var.cxx) - part E: the
//! coefficient engines (mma2cfu_, mma2cfv_, mma2er1_, mma2er2_) and the
//! approximation drivers (mma2ce1_, mma2ce2_).
//!
//! Encoding note: as in parts A-D, the OCCT Fortran 1-based index formulas
//! are kept on 0-based slices by subtracting the OCCT virtual origin inside
//! the flat-index expression; the ce1_ workspace zones are expressed with
//! split_at_mut over the raw allocation.
//!
//! OCCT quirk preserved (AdvApp2Var_ApproxF2var.cxx L3952-3953): both
//! mma2jmx_ results are written to the XMAXJV zone (&wrkar_off[ipt5]); the
//! XMAXJU zone (&wrkar_off[ipt4]) is never initialized by mma2ce1_ and is
//! read downstream.  OCCT reads uninitialized heap there; the rcad encoding
//! keeps deterministic zeros (fresh Vec allocations are zeroed).

use super::approxf2var_b::mma2jmx_;
use super::approxf2var_b::mmapptt_;
use super::approxf2var_b::mma2moy_;
use super::math_base_b::mzsnorm_;
use super::sys_base;

// ---------------------------------------------------------------------------
// mma2cfu_ (AdvApp2Var_ApproxF2var.cxx L4998-5226)
// ---------------------------------------------------------------------------

/// OCCT mma2cfu_ (AdvApp2Var_ApproxF2var.cxx L4998-5226) - terms connected to
/// degree NDUJAC by U.  sosotb / diditb / gssutb / chpair are raw 0-based
/// sub-buffers; soditb / disotb are adjusted sub-buffers (offset dim1 + 1);
/// chimpr is 1-based emulated.
#[allow(clippy::too_many_arguments)]
pub fn mma2cfu_(
    ndujac: &mut i32,
    nbpntu: &i32,
    nbpntv: &i32,
    sosotb: &[f64],
    disotb: &[f64],
    soditb: &[f64],
    diditb: &[f64],
    gssutb: &[f64],
    chpair: &mut [f64],
    chimpr: &mut [f64],
) {
    // OCCT: soditb_dim1 = *nbpntu / 2; soditb_offset = dim1 + 1.
    let soditb_dim1 = (nbpntu / 2) as i32;
    let disotb_dim1 = (nbpntu / 2) as i32;
    // OCCT: sosotb_dim1 = *nbpntu / 2 + 1 (raw, no adjustment).
    let sosotb_dim1 = (nbpntu / 2 + 1) as i32;
    let diditb_dim1 = (nbpntu / 2 + 1) as i32;

    let sd = |ii: i32, jj: i32| -> usize {
        // OCCT: soditb[ii + jj * dim1] - (dim1 + 1).
        ((ii - 1) + (jj - 1) * soditb_dim1) as usize
    };
    let di = |ii: i32, jj: i32| -> usize {
        ((ii - 1) + (jj - 1) * disotb_dim1) as usize
    };
    let ss = |ii: i32, jj: i32| -> usize {
        // OCCT: sosotb[ii + jj * dim1] (raw, 0-based).
        (ii + jj * sosotb_dim1) as usize
    };
    let dd = |ii: i32, jj: i32| -> usize {
        // OCCT: diditb[ii + jj * dim1] (raw, 0-based).
        (ii + jj * diditb_dim1) as usize
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2CFU");
    }

    let nptu2 = *nbpntu / 2;
    let nptv2 = *nbpntv / 2;

    // ----------------- Calculate  coefficients of even degree --------------
    if *ndujac % 2 == 0 {
        for jj in 1..=nptv2 {
            let mut bid1 = 0.;
            let mut bid2 = 0.;
            for ii in 1..=nptu2 {
                let bid0 = gssutb[ii as usize];
                bid1 += sosotb[ss(ii, jj)] * bid0;
                bid2 += soditb[sd(ii, jj)] * bid0;
            }
            // OCCT: chpair[jj] (raw 0-based); chimpr[jj] (1-based emulated).
            chpair[jj as usize] = bid1;
            chimpr[(jj - 1) as usize] = bid2;
        }
    // --------------- Calculate coefficients of uneven degree ----------
    } else {
        for jj in 1..=nptv2 {
            let mut bid1 = 0.;
            let mut bid2 = 0.;
            for ii in 1..=nptu2 {
                let bid0 = gssutb[ii as usize];
                bid1 += disotb[di(ii, jj)] * bid0;
                bid2 += diditb[dd(ii, jj)] * bid0;
            }
            chpair[jj as usize] = bid1;
            chimpr[(jj - 1) as usize] = bid2;
        }
    }

    // ------- Add terms connected to the supplementary root (0.D0) ------ */
    // ----------- of Legendre polynom of uneven degree NBPNTU -----------
    // --> Only even NDUJAC terms are modified as GSSUTB(0) = 0
    //     when NDUJAC is uneven.
    if *nbpntu % 2 != 0 && *ndujac % 2 == 0 {
        let bid0 = gssutb[0];
        for jj in 1..=nptv2 {
            chpair[jj as usize] += sosotb[ss(0, jj)] * bid0;
            chimpr[(jj - 1) as usize] += diditb[dd(0, jj)] * bid0;
        }
    }

    // ------ Calculate the terms connected to supplementary roots (0.D0) ------
    // ----------- of Legendre polynom of uneven degree NBPNTV -----------
    if *nbpntv % 2 != 0 {
        // --> Only CHPAIR terms are calculated as GSSVTB(0,IH-IDEBV)=0
        //     when IH is uneven (see MMA2CFV).
        if *ndujac % 2 == 0 {
            let mut bid1 = 0.;
            for ii in 1..=nptu2 {
                bid1 += sosotb[ss(ii, 0)] * gssutb[ii as usize];
            }
            chpair[0] = bid1;
        } else {
            let mut bid1 = 0.;
            for ii in 1..=nptu2 {
                bid1 += diditb[dd(ii, 0)] * gssutb[ii as usize];
            }
            chpair[0] = bid1;
        }
        if *nbpntu % 2 != 0 {
            chpair[0] += sosotb[0] * gssutb[0];
        }
    }

    // ------------------------------ The end -------------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2CFU");
    }
}

// ---------------------------------------------------------------------------
// mma2cfv_ (AdvApp2Var_ApproxF2var.cxx L5228-5381)
// ---------------------------------------------------------------------------

/// OCCT mma2cfv_ (AdvApp2Var_ApproxF2var.cxx L5228-5381) - coefficients of
/// degree NDVJAC by V, degrees MINDGU..MAXDGU by U.  gssvtb is a raw
/// 0-based sub-buffer; chpair / chimpr are adjusted sub-buffers relative to
/// MINDGU; patjac is the adjusted 1-based emulated output sub-buffer.
#[allow(clippy::too_many_arguments)]
pub fn mma2cfv_(
    ndvjac: &mut i32,
    mindgu: &mut i32,
    maxdgu: &mut i32,
    nbpntv: &i32,
    gssvtb: &[f64],
    chpair: &[f64],
    chimpr: &[f64],
    patjac: &mut [f64],
) {
    // OCCT: patjac_offset = *mindgu.
    let patjac_off = *mindgu;
    // OCCT: chimpr_dim1 = *nbpntv / 2; chimpr_offset = dim1 * mindgu + 1.
    let chimpr_dim1 = (nbpntv / 2) as i32;
    let chimpr_off = chimpr_dim1 * *mindgu + 1;
    // OCCT: chpair_dim1 = *nbpntv / 2 + 1; chpair_offset = dim1 * mindgu.
    let chpair_dim1 = (nbpntv / 2 + 1) as i32;
    let chpair_off = chpair_dim1 * *mindgu;

    let cp = |jj: i32, ii: i32| -> usize {
        // OCCT: chpair[jj + ii * dim1] - dim1 * mindgu.
        ((jj - 1) + (ii - *mindgu) * chpair_dim1) as usize
    };
    let cm = |jj: i32, ii: i32| -> usize {
        // OCCT: chimpr[jj + ii * dim1] - (dim1 * mindgu + 1).
        ((jj - 1) + (ii - *mindgu) * chimpr_dim1) as usize
    };
    let pj = |ii: i32| -> usize {
        // OCCT: patjac[ii] - mindgu.
        (ii - patjac_off) as usize
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2CFV");
    }
    let nptv2 = *nbpntv / 2;

    // --------- Calculate the coefficients for even degree NDVJAC ----------
    if *ndvjac % 2 == 0 {
        for ii in *mindgu..=*maxdgu {
            let mut bid1 = 0.;
            for jj in 1..=nptv2 {
                bid1 += chpair[cp(jj, ii)] * gssvtb[jj as usize];
            }
            patjac[pj(ii)] = bid1;
        }
    // -------- Calculate the coefficients for uneven degree NDVJAC -----
    } else {
        for ii in *mindgu..=*maxdgu {
            let mut bid1 = 0.;
            for jj in 1..=nptv2 {
                bid1 += chimpr[cm(jj, ii)] * gssvtb[jj as usize];
            }
            patjac[pj(ii)] = bid1;
        }
    }

    // ------- Add terms connected to the supplementary root (0.D0) ----- */
    // --------of the Legendre polynom of uneven degree  NBPNTV ---------
    if *nbpntv % 2 != 0 && *ndvjac % 2 == 0 {
        let bid1 = gssvtb[0];
        for ii in *mindgu..=*maxdgu {
            // OCCT: patjac[ii] += bid1 * chpair[ii * chpair_dim1].
            patjac[pj(ii)] += bid1 * chpair[((ii - *mindgu) * chpair_dim1) as usize];
        }
    }

    // ------------------------------ The end -------------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2CFV");
    }
}

// ---------------------------------------------------------------------------
// mma2er1_ (AdvApp2Var_ApproxF2var.cxx L6209-6372)
// ---------------------------------------------------------------------------

/// OCCT mma2er1_ (AdvApp2Var_ApproxF2var.cxx L6209-6372) - max approximation
/// error done when the PATJAC coefficients in [MINDGU, MAXDGU] x
/// [MINDGV, MAXDGV] are removed.  patjac / xmaxju / xmaxjv / vecerr are raw
/// sub-buffers (vecerr 1-based emulated).
#[allow(clippy::too_many_arguments)]
pub fn mma2er1_(
    ndjacu: &mut i32,
    ndjacv: &mut i32,
    ndimen: &mut i32,
    mindgu: &mut i32,
    maxdgu: &mut i32,
    mindgv: &mut i32,
    maxdgv: &mut i32,
    iordru: &mut i32,
    iordrv: &mut i32,
    xmaxju: &[f64],
    xmaxjv: &[f64],
    patjac: &[f64],
    vecerr: &mut [f64],
    erreur: &mut f64,
) {
    // OCCT: patjac_dim1 = *ndjacu + 1; patjac_dim2 = *ndjacv + 1;
    //       patjac_offset = dim1 * dim2.
    let patjac_dim1 = (*ndjacu + 1) as usize;
    let patjac_dim2 = (*ndjacv + 1) as usize;

    let pj = |ii: i32, jj: i32, nd: i32| -> usize {
        // OCCT: patjac[ii + (jj + nd * dim2) * dim1] - dim1 * dim2.
        ((ii - 1) + ((jj - 1) + (nd - 1) * patjac_dim2 as i32) * patjac_dim1 as i32) as usize
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2ER1");
    }

    let minu = (*iordru + 1) << 1;
    let minv = (*iordrv + 1) << 1;

    // ------------------- Calculate the increment of the max error --------
    for nd in 1..=*ndimen {
        let mut bid1 = 0.;
        for jj in *mindgv..=*maxdgv {
            let mut bid0 = 0.;
            for ii in *mindgu..=*maxdgu {
                bid0 += patjac[pj(ii, jj, nd)].abs() * xmaxju[(ii - minu) as usize];
            }
            bid1 = bid0 * xmaxjv[(jj - minv) as usize] + bid1;
        }
        // OCCT: vecerr[nd] (1-based emulated).
        vecerr[(nd - 1) as usize] = bid1;
    }

    // ----------------------- Calculate the max error  --------------------
    let bid1 = mzsnorm_(ndimen, vecerr);
    let mut vaux = [0.0f64; 2];
    vaux[0] = *erreur;
    vaux[1] = bid1;
    let mut nd = 2;
    *erreur = mzsnorm_(&mut nd, &mut vaux);

    if ldbg {
        sys_base::mgsomsg_("MMA2ER1");
    }
}

// ---------------------------------------------------------------------------
// mma2er2_ (AdvApp2Var_ApproxF2var.cxx L6374-6622)
// ---------------------------------------------------------------------------

/// OCCT mma2er2_ (AdvApp2Var_ApproxF2var.cxx L6374-6622) - remove PATJAC
/// coefficients to reach the minimum degree by U and V checking the
/// tolerance.  Buffer conventions as mma2er1_.
#[allow(clippy::too_many_arguments)]
pub fn mma2er2_(
    ndjacu: &mut i32,
    ndjacv: &mut i32,
    ndimen: &mut i32,
    mindgu: &mut i32,
    maxdgu: &mut i32,
    mindgv: &mut i32,
    maxdgv: &mut i32,
    iordru: &mut i32,
    iordrv: &mut i32,
    xmaxju: &[f64],
    xmaxjv: &[f64],
    patjac: &[f64],
    epmscut: &mut f64,
    vecerr: &mut [f64],
    erreur: &mut f64,
    newdgu: &mut i32,
    newdgv: &mut i32,
) {
    let patjac_dim1 = (*ndjacu + 1) as usize;
    let patjac_dim2 = (*ndjacv + 1) as usize;

    let pj = |ii: i32, jj: i32, nd: i32| -> usize {
        ((ii - 1) + ((jj - 1) + (nd - 1) * patjac_dim2 as i32) * patjac_dim1 as i32) as usize
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2ER2");
    }

    let i2rdu = (*iordru + 1) << 1;
    let i2rdv = (*iordrv + 1) << 1;
    let mut nu = *maxdgu;
    let mut nv = *maxdgv;

    // -------------------- Cutting of coefficients ------------------------
    // OCCT L1001 loop:
    loop {
        // ------------------- Calculate the increment of max error --------
        // ----- during the removal of coeff. of indices MINDGU..MAXDGU ----
        // ---------------- by U, the degree by V is fixed to NV -----------
        if nv > *mindgv {
            let bid0 = xmaxjv[(nv - i2rdv) as usize];
            for nd in 1..=*ndimen {
                let mut bid1 = 0.;
                for ii in i2rdu..=nu {
                    bid1 += patjac[pj(ii, nv, nd)].abs()
                        * xmaxju[(ii - i2rdu) as usize]
                        * bid0;
                }
                vecerr[(nd - 1) as usize] = bid1;
            }
        } else {
            vecerr[0] = *epmscut * 2.;
        }
        let mut errnv = mzsnorm_(ndimen, vecerr);

        // ------------------- Calculate the increment of max error --------
        // ----- during the removal of coeff. of indices MINDGV..MAXDGV ----
        // ---------------- by V, the degree by U is fixed to NU -----------
        if nu > *mindgu {
            let bid0 = xmaxju[(nu - i2rdu) as usize];
            for nd in 1..=*ndimen {
                let mut bid1 = 0.;
                for jj in i2rdv..=nv {
                    bid1 += patjac[pj(nu, jj, nd)].abs()
                        * xmaxjv[(jj - i2rdv) as usize]
                        * bid0;
                }
                vecerr[(nd - 1) as usize] = bid1;
            }
        } else {
            vecerr[0] = *epmscut * 2.;
        }
        let mut errnu = mzsnorm_(ndimen, vecerr);

        // ----------------------- Calculate the max error -----------------
        let mut vaux = [0.0f64; 2];
        vaux[0] = *erreur;
        vaux[1] = errnu;
        let mut nd = 2;
        errnu = mzsnorm_(&mut nd, &mut vaux);
        vaux[1] = errnv;
        errnv = mzsnorm_(&mut nd, &mut vaux);

        if errnu > errnv {
            if errnv < *epmscut {
                *erreur = errnv;
                nv -= 1;
            } else {
                break; // OCCT goto L2001.
            }
        } else {
            if errnu < *epmscut {
                *erreur = errnu;
                nu -= 1;
            } else {
                break; // OCCT goto L2001.
            }
        }
    }

    // -------------------------- Return the degrees -----------------------
    // OCCT L2001:
    *newdgu = nu.max(1);
    *newdgv = nv.max(1);

    // ----------------------------------- The end -------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2ER2");
    }
}

// ---------------------------------------------------------------------------
// mma2ce2_ (AdvApp2Var_ApproxF2var.cxx L4023-4996)
// ---------------------------------------------------------------------------

/// OCCT mma2ce2_ (AdvApp2Var_ApproxF2var.cxx L4023-4996) - calculation of
/// the approximation coefficients by zones and of the cutting tests.  All
/// buffers are raw sub-buffers handed by mma2ce1_.
#[allow(clippy::too_many_arguments)]
pub fn mma2ce2_(
    numdec: &mut i32,
    ndimen: &mut i32,
    nbsesp: &mut i32,
    ndimse: &[i32],
    ndminu: &mut i32,
    ndminv: &mut i32,
    ndguli: &mut i32,
    ndgvli: &mut i32,
    ndjacu: &mut i32,
    ndjacv: &mut i32,
    iordru: &mut i32,
    iordrv: &mut i32,
    nbpntu: &mut i32,
    nbpntv: &mut i32,
    epsapr: &[f64],
    sosotb: &mut [f64],
    disotb: &mut [f64],
    soditb: &mut [f64],
    diditb: &mut [f64],
    gssutb: &mut [f64],
    gssvtb: &mut [f64],
    xmaxju: &mut [f64],
    xmaxjv: &mut [f64],
    vecerr: &mut [f64],
    chpair: &mut [f64],
    chimpr: &mut [f64],
    patjac: &mut [f64],
    errmax: &mut [f64],
    errmoy: &mut [f64],
    ndegpu: &mut i32,
    ndegpv: &mut i32,
    itydec: &mut i32,
    iercod: &mut i32,
) {
    // OCCT: vecerr_dim1 = *ndimen (1-based emulated buffer).
    let vecerr_dim1 = *ndimen;
    // OCCT: patjac_dim1 = *ndjacu + 1; patjac_dim2 = *ndjacv + 1.
    let patjac_dim1 = (*ndjacu + 1) as i32;
    let patjac_dim2 = (*ndjacv + 1) as i32;
    // OCCT: gssutb_dim1 = *nbpntu / 2 + 1; gssvtb_dim1 = *nbpntv / 2 + 1.
    let gssutb_dim1 = (*nbpntu / 2 + 1) as i32;
    let gssvtb_dim1 = (*nbpntv / 2 + 1) as i32;
    // OCCT: chimpr_dim1 = *nbpntv / 2; chimpr_dim2 = *ndjacu - 2*(iordru+1) + 1.
    let chimpr_dim1 = (*nbpntv / 2) as i32;
    let chimpr_dim2 = *ndjacu - ((*iordru + 1) << 1) + 1;
    // OCCT: chpair_dim1 = *nbpntv / 2 + 1; chpair_dim2 = chimpr_dim2.
    let chpair_dim1 = (*nbpntv / 2 + 1) as i32;
    let chpair_dim2 = chimpr_dim2;
    // sosotb family dims (adjusted buffers).
    let sosotb_dim1 = (*nbpntu / 2 + 1) as i32;
    let sosotb_dim2 = (*nbpntv / 2 + 1) as i32;
    let disotb_dim1 = (*nbpntu / 2) as i32;
    let disotb_dim2 = (*nbpntv / 2) as i32;
    let soditb_dim1 = (*nbpntu / 2) as i32;
    let soditb_dim2 = (*nbpntv / 2) as i32;
    let diditb_dim1 = (*nbpntu / 2 + 1) as i32;
    let diditb_dim2 = (*nbpntv / 2 + 1) as i32;

    let ve = |slot: i32, nd: i32| -> usize {
        // OCCT: vecerr[nd + slot * vecerr_dim1] - (dim1 + 1).
        ((nd - 1) + (slot - 1) * vecerr_dim1) as usize
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2CE2");
    }
    // --> A priori everything is OK
    *iercod = 0;
    // --> test of inputs
    if *numdec < 0 || *numdec > 5 || ((*iordru << 1) + 1 > *ndminu) || *ndminu > *ndguli
        || *ndguli >= *ndjacu || ((*iordrv << 1) + 1 > *ndminv) || *ndminv > *ndgvli
        || *ndgvli >= *ndjacv
    {
        // OCCT goto L9001.
        *iercod = 1;
        sys_base::maermsg_("MMA2CE2", iercod);
        if ldbg {
            sys_base::mgsomsg_("MMA2CE2");
        }
        return;
    }
    // --> A priori, no cuts to be done
    *itydec = 0;
    // --> Min. degrees to return: NDMINU,NDMINV
    *ndegpu = *ndminu;
    *ndegpv = *ndminv;
    // --> For the moment, max errors are null
    sys_base::mvriraz_(*nbsesp, errmax);
    let mut nd = *ndimen << 2;
    sys_base::mvriraz_(nd, vecerr);
    // --> and the square, too.
    nd = (*ndjacu + 1) * (*ndjacv + 1) * *ndimen;
    sys_base::mvriraz_(nd, patjac);

    let i2rdu = (*iordru + 1) << 1;
    let i2rdv = (*iordrv + 1) << 1;

    // Misc locals (the goto structure is reproduced with labelled blocks).
    let mut minu;
    let mut maxu;
    let mut minv;
    let mut maxv;
    let mut igsu;
    let mut igsv;
    let mut idim;
    let mut ndses;
    let mut nu = 0;
    let mut nv = 0;
    let mut vaux = [0.0f64; 3];
    let mut ii3 = 0;

    // **********************************************************************
    // -------------------- HERE IT IS POSSIBLE TO CUT ----------------------
    // **********************************************************************
    if *numdec > 0 && *numdec <= 5 {
        // ---------------------- Calculate coeff of zone 4 ----------------
        minu = *ndguli + 1;
        maxu = *ndjacu;
        minv = *ndgvli + 1;
        maxv = *ndjacv;
        if minu > maxu || minv > maxv {
            // OCCT goto L9001.
            *iercod = 1;
            sys_base::maermsg_("MMA2CE2", iercod);
            if ldbg {
                sys_base::mgsomsg_("MMA2CE2");
            }
            return;
        }

        // ---------------- Calculate the terms connected to degree by U ---
        for nd in 1..=*ndimen {
            for kk in minu..=maxu {
                igsu = kk - i2rdu;
                // OCCT: mma2cfu_(&kk, nbpntu, nbpntv,
                //   &sosotb[nd * dim2 * dim1], &disotb[(nd*dim2+1)*dim1+1],
                //   &soditb[(nd*dim2+1)*dim1+1], &diditb[nd*dim2*dim1],
                //   &gssutb[igsu*dim1], &chpair[(igsu+nd*dim2)*dim1],
                //   &chimpr[(igsu+nd*dim2)*dim1+1]) (all sub-pointers raw).
                let mut kk_ = kk;
                mma2cfu_(
                    &mut kk_,
                    nbpntu,
                    nbpntv,
                    &sosotb[((nd - 1) * sosotb_dim2 * sosotb_dim1) as usize..],
                    &disotb[((nd - 1) * disotb_dim2 * disotb_dim1) as usize..],
                    &soditb[((nd - 1) * soditb_dim2 * soditb_dim1) as usize..],
                    &diditb[((nd - 1) * diditb_dim2 * diditb_dim1) as usize..],
                    &gssutb[(igsu * gssutb_dim1) as usize..],
                    &mut chpair[((igsu + (nd - 1) * chpair_dim2) * chpair_dim1) as usize..],
                    &mut chimpr[
                        ((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..
                    ],
                );
            }
        }

        // ------------------- Calculate the coefficients of PATJAC --------
        igsu = minu - i2rdu;
        let mut jj = minv;
        while jj <= maxv {
            igsv = jj - i2rdv;
            for nd in 1..=*ndimen {
                let mut jj_ = jj;
                let mut minu_ = minu;
                let mut maxu_ = maxu;
                mma2cfv_(
                    &mut jj_,
                    &mut minu_,
                    &mut maxu_,
                    nbpntv,
                    &gssvtb[(igsv * gssvtb_dim1) as usize..],
                    &chpair[((igsu + (nd - 1) * chpair_dim2) * chpair_dim1) as usize..],
                    &chimpr[((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..],
                    &mut patjac[((minu - 1)
                        + ((jj - 1) + (nd - 1) * patjac_dim2) * patjac_dim1)
                        as usize..],
                );
            }

            // ----- Contribution of calculated terms to the error ---------
            // for terms (I,J) with MINU <= I <= MAXU, J fixe.
            idim = 1;
            for nd in 1..=*nbsesp {
                ndses = ndimse[(nd - 1) as usize];
                let mut ndses_ = ndses;
                let mut minu_ = minu;
                let mut maxu_ = maxu;
                let mut jj_ = jj;
                let mut jj2_ = jj;
                let mut err = 0.;
                mma2er1_(
                    ndjacu,
                    ndjacv,
                    &mut ndses_,
                    &mut minu_,
                    &mut maxu_,
                    &mut jj_,
                    &mut jj2_,
                    iordru,
                    iordrv,
                    xmaxju,
                    xmaxjv,
                    &patjac
                        [((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                    &mut vecerr[..],
                    &mut err,
                );
                vecerr[ve(4, nd)] = err;
                if vecerr[ve(4, nd)] > epsapr[(nd - 1) as usize] {
                    // OCCT goto L9300.
                    itydec_from_numdec(numdec, itydec, iercod, 9300, ldbg);
                    return;
                }
                idim += ndses;
            }
            jj += 1;
        }

        // ---------------------- Calculate the coeff of zone 2 ------------
        minu = (*iordru + 1) << 1;
        maxu = *ndguli;
        minv = *ndgvli + 1;
        maxv = *ndjacv;

        // --> If zone 2 is empty, pass to zone 3.
        //     VECERR(ND,2) was already set to zero.
        if minu <= maxu {
            // ---------------- Calculate the terms connected to degree by U
            for nd in 1..=*ndimen {
                for kk in minu..=maxu {
                    igsu = kk - i2rdu;
                    let mut kk_ = kk;
                    mma2cfu_(
                        &mut kk_,
                        nbpntu,
                        nbpntv,
                        &sosotb[((nd - 1) * sosotb_dim2 * sosotb_dim1) as usize..],
                        &disotb[((nd - 1) * disotb_dim2 * disotb_dim1) as usize..],
                        &soditb[((nd - 1) * soditb_dim2 * soditb_dim1) as usize..],
                        &diditb[((nd - 1) * diditb_dim2 * diditb_dim1) as usize..],
                        &gssutb[(igsu * gssutb_dim1) as usize..],
                        &mut chpair
                            [((igsu + (nd - 1) * chpair_dim2) * chpair_dim1) as usize..],
                        &mut chimpr[
                            ((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..
                        ],
                    );
                }
            }

            // ------------------- Calculate the coefficients of PATJAC ----
            igsu = minu - i2rdu;
            for jj in minv..=maxv {
                igsv = jj - i2rdv;
                for nd in 1..=*ndimen {
                    let mut jj_ = jj;
                    let mut minu_ = minu;
                    let mut maxu_ = maxu;
                    mma2cfv_(
                        &mut jj_,
                        &mut minu_,
                        &mut maxu_,
                        nbpntv,
                        &gssvtb[(igsv * gssvtb_dim1) as usize..],
                        &chpair
                            [((igsu + (nd - 1) * chpair_dim2) * chpair_dim1) as usize..],
                        &chimpr[
                            ((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..
                        ],
                        &mut patjac[((minu - 1)
                            + ((jj - 1) + (nd - 1) * patjac_dim2) * patjac_dim1)
                            as usize..],
                    );
                }
            }

            // -----Contribution of calculated terms to the error ---------
            // for terms (I,J) with MINU <= I <= MAXU, MINV <= J <= MAXV
            idim = 1;
            for nd in 1..=*nbsesp {
                ndses = ndimse[(nd - 1) as usize];
                let mut ndses_ = ndses;
                let mut minu_ = minu;
                let mut maxu_ = maxu;
                let mut minv_ = minv;
                let mut maxv_ = maxv;
                let mut err = 0.;
                mma2er1_(
                    ndjacu,
                    ndjacv,
                    &mut ndses_,
                    &mut minu_,
                    &mut maxu_,
                    &mut minv_,
                    &mut maxv_,
                    iordru,
                    iordrv,
                    xmaxju,
                    xmaxjv,
                    &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                    &mut vecerr[..],
                    &mut err,
                );
                vecerr[ve(2, nd)] = err;
                idim += ndses;
            }
        }

        // ---------------------- Calculation of coeff of zone 3 -----------
        // OCCT L300:
        minu = *ndguli + 1;
        maxu = *ndjacu;
        minv = (*iordrv + 1) << 1;
        maxv = *ndgvli;

        // -> If zone 3 is empty, pass to the test of cutting.
        //    VECERR(ND,3) was already set to zero
        if minv <= maxv {
            // ---- The terms connected to degree by U are already calculated
            // ------------------- Calculation of coefficients of PATJAC ----
            igsu = minu - i2rdu;
            for jj in minv..=maxv {
                igsv = jj - i2rdv;
                for nd in 1..=*ndimen {
                    let mut jj_ = jj;
                    let mut minu_ = minu;
                    let mut maxu_ = maxu;
                    mma2cfv_(
                        &mut jj_,
                        &mut minu_,
                        &mut maxu_,
                        nbpntv,
                        &gssvtb[(igsv * gssvtb_dim1) as usize..],
                        &chpair[((igsu + (nd - 1) * chpair_dim2) * chpair_dim1)
                            as usize..],
                        &chimpr[
                            ((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..
                        ],
                        &mut patjac[((minu - 1)
                            + ((jj - 1) + (nd - 1) * patjac_dim2) * patjac_dim1)
                            as usize..],
                    );
                }
            }

            // ----- Contribution of calculated terms to the error ---------
            // for terms (I,J) with MINU <= I <= MAXU, MINV <= J <= MAXV.
            idim = 1;
            for nd in 1..=*nbsesp {
                ndses = ndimse[(nd - 1) as usize];
                let mut ndses_ = ndses;
                let mut minu_ = minu;
                let mut maxu_ = maxu;
                let mut minv_ = minv;
                let mut maxv_ = maxv;
                let mut err = 0.;
                mma2er1_(
                    ndjacu,
                    ndjacv,
                    &mut ndses_,
                    &mut minu_,
                    &mut maxu_,
                    &mut minv_,
                    &mut maxv_,
                    iordru,
                    iordrv,
                    xmaxju,
                    xmaxjv,
                    &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                    &mut vecerr[..],
                    &mut err,
                );
                vecerr[ve(3, nd)] = err;
                idim += ndses;
            }
        }

        // --------------------------- Tests of cutting --------------------
        // OCCT L400:
        for nd in 1..=*nbsesp {
            vaux[0] = vecerr[ve(2, nd)];
            vaux[1] = vecerr[ve(4, nd)];
            vaux[2] = vecerr[ve(3, nd)];
            ii3 = 3;
            let mut vaux3 = vaux;
            errmax[(nd - 1) as usize] = mzsnorm_(&mut ii3, &mut vaux3);
            if errmax[(nd - 1) as usize] > epsapr[(nd - 1) as usize] {
                let mut ii2 = 2;
                let mut zv_vaux = vaux;
                let zv = mzsnorm_(&mut ii2, &mut zv_vaux);
                let zu_vaux = [vaux[1], vaux[2]];
                let mut zu_vaux = zu_vaux;
                let zu = mzsnorm_(&mut ii2, &mut zu_vaux);
                if zu > epsapr[(nd - 1) as usize] && zv > epsapr[(nd - 1) as usize] {
                    // OCCT goto L9300.
                    itydec_from_numdec(numdec, itydec, iercod, 9300, ldbg);
                    return;
                }
                if zu > zv {
                    // OCCT goto L9100.
                    itydec_from_numdec(numdec, itydec, iercod, 9100, ldbg);
                    return;
                } else {
                    // OCCT goto L9200.
                    itydec_from_numdec(numdec, itydec, iercod, 9200, ldbg);
                    return;
                }
            }
        }

        // --- OK, the square is valid, the coeff of zone 1 are calculated -
        minu = (*iordru + 1) << 1;
        maxu = *ndguli;
        minv = (*iordrv + 1) << 1;
        maxv = *ndgvli;

        // --> If zone 1 is empty, pass to the calculation of Max and Average
        if minu <= maxu && minv <= maxv {
            // ---- The terms connected to degree by U are already calculated
            // ------------------- Calculate the coefficients of PATJAC ----
            igsu = minu - i2rdu;
            for jj in minv..=maxv {
                igsv = jj - i2rdv;
                for nd in 1..=*ndimen {
                    let mut jj_ = jj;
                    let mut minu_ = minu;
                    let mut maxu_ = maxu;
                    mma2cfv_(
                        &mut jj_,
                        &mut minu_,
                        &mut maxu_,
                        nbpntv,
                        &gssvtb[(igsv * gssvtb_dim1) as usize..],
                        &chpair[((igsu + (nd - 1) * chpair_dim2) * chpair_dim1)
                            as usize..],
                        &chimpr[
                            ((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..
                        ],
                        &mut patjac[((minu - 1)
                            + ((jj - 1) + (nd - 1) * patjac_dim2) * patjac_dim1)
                            as usize..],
                    );
                }
            }
        }

        // --------------- Now the degree is maximally lowered ------------
        // OCCT L600:
        minu = 1.max(((*iordru << 1) + 1)).max(*ndminu);
        maxu = *ndguli;
        minv = 1.max(((*iordrv << 1) + 1)).max(*ndminv);
        maxv = *ndgvli;
        idim = 1;
        for nd in 1..=*nbsesp {
            ndses = ndimse[(nd - 1) as usize];
            if maxu >= (*iordru + 1) << 1 && maxv >= (*iordrv + 1) << 1 {
                let mut ndses_ = ndses;
                let mut minu_ = minu;
                let mut maxu_ = maxu;
                let mut minv_ = minv;
                let mut maxv_ = maxv;
                let mut epsapr_nd = epsapr[(nd - 1) as usize];
                let mut err_nd = errmax[(nd - 1) as usize];
                mma2er2_(
                    ndjacu,
                    ndjacv,
                    &mut ndses_,
                    &mut minu_,
                    &mut maxu_,
                    &mut minv_,
                    &mut maxv_,
                    iordru,
                    iordrv,
                    xmaxju,
                    xmaxjv,
                    &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                    &mut epsapr_nd,
                    &mut vecerr[..],
                    &mut err_nd,
                    &mut nu,
                    &mut nv,
                );
                errmax[(nd - 1) as usize] = err_nd;
            } else {
                nu = maxu;
                nv = maxv;
            }
            let nu1 = nu + 1;
            let nv1 = nv + 1;

            // --> Calculate the average error.
            let mut ndses_ = ndses;
            let mut nu1_ = nu1;
            let mut nv1_ = nv1;
            let mut ndjacu_a = *ndjacu;
            let mut ndjacu_b = *ndjacu;
            let mut ndjacv_a = *ndjacv;
            let mut ndjacv_b = *ndjacv;
            let mut errmoy_nd = errmoy[(nd - 1) as usize];
            mma2moy_(
                &mut ndjacu_a,
                &mut ndjacv_a,
                &mut ndses_,
                &mut nu1_,
                &mut ndjacu_b,
                &mut nv1_,
                &mut ndjacv_b,
                iordru,
                iordrv,
                &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                &mut errmoy_nd,
            );
            errmoy[(nd - 1) as usize] = errmoy_nd;

            // --> Set to 0.D0 the rejected coeffs.
            for ii in idim..=idim + ndses - 1 {
                for jj in nv1..=*ndjacv {
                    for kk in nu1..=*ndjacu {
                        patjac[((kk - 1)
                            + ((jj - 1) + (ii - 1) * patjac_dim2) * patjac_dim1)
                            as usize] = 0.;
                    }
                }
            }

            // --> Return the nb of coeffs of approximation.
            *ndegpu = (*ndegpu).max(nu);
            *ndegpv = (*ndegpv).max(nv);
            idim += ndses;
        }
    } else {
        // -------------------- IT IS NOT POSSIBLE TO CUT ------------------
        minu = (*iordru + 1) << 1;
        maxu = *ndjacu;
        minv = (*iordrv + 1) << 1;
        maxv = *ndjacv;

        // ---------------- Calculate the terms connected to degree by U ---
        for nd in 1..=*ndimen {
            for kk in minu..=maxu {
                igsu = kk - i2rdu;
                let mut kk_ = kk;
                mma2cfu_(
                    &mut kk_,
                    nbpntu,
                    nbpntv,
                    &sosotb[((nd - 1) * sosotb_dim2 * sosotb_dim1) as usize..],
                    &disotb[((nd - 1) * disotb_dim2 * disotb_dim1) as usize..],
                    &soditb[((nd - 1) * soditb_dim2 * soditb_dim1) as usize..],
                    &diditb[((nd - 1) * diditb_dim2 * diditb_dim1) as usize..],
                    &gssutb[(igsu * gssutb_dim1) as usize..],
                    &mut chpair[((igsu + (nd - 1) * chpair_dim2) * chpair_dim1) as usize..],
                    &mut chimpr[((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..],
                );
            }

            // ---------------------- Calculate all coefficients ------------
            igsu = minu - i2rdu;
            for jj in minv..=maxv {
                igsv = jj - i2rdv;
                let mut jj_ = jj;
                let mut minu_ = minu;
                let mut maxu_ = maxu;
                mma2cfv_(
                    &mut jj_,
                    &mut minu_,
                    &mut maxu_,
                    nbpntv,
                    &gssvtb[(igsv * gssvtb_dim1) as usize..],
                    &chpair[((igsu + (nd - 1) * chpair_dim2) * chpair_dim1) as usize..],
                    &chimpr[((igsu + (nd - 1) * chimpr_dim2) * chimpr_dim1) as usize..],
                    &mut patjac[((minu - 1)
                        + ((jj - 1) + (nd - 1) * patjac_dim2) * patjac_dim1)
                        as usize..],
                );
            }
        }

        // ----- Contribution of calculated terms to the approximation error
        // for terms (I,J) with MINU <= I <= MAXU, MINV <= J <= MAXV
        idim = 1;
        for nd in 1..=*nbsesp {
            ndses = ndimse[(nd - 1) as usize];
            minu = (*iordru + 1) << 1;
            maxu = *ndjacu;
            minv = *ndgvli + 1;
            maxv = *ndjacv;
            let mut ndses_ = ndses;
            let mut minu_ = minu;
            let mut maxu_ = maxu;
            let mut minv_ = minv;
            let mut maxv_ = maxv;
            let mut err = 0.;
            mma2er1_(
                ndjacu,
                ndjacv,
                &mut ndses_,
                &mut minu_,
                &mut maxu_,
                &mut minv_,
                &mut maxv_,
                iordru,
                iordrv,
                xmaxju,
                xmaxjv,
                &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                &mut vecerr[..],
                &mut err,
            );
            errmax[(nd - 1) as usize] = err;
            minu = *ndguli + 1;
            maxu = *ndjacu;
            minv = (*iordrv + 1) << 1;
            maxv = *ndgvli;
            if minv <= maxv {
                let mut ndses_ = ndses;
                let mut minu_ = minu;
                let mut maxu_ = maxu;
                let mut minv_ = minv;
                let mut maxv_ = maxv;
                let mut err = errmax[(nd - 1) as usize];
                mma2er1_(
                    ndjacu,
                    ndjacv,
                    &mut ndses_,
                    &mut minu_,
                    &mut maxu_,
                    &mut minv_,
                    &mut maxv_,
                    iordru,
                    iordrv,
                    xmaxju,
                    xmaxjv,
                    &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                    &mut vecerr[..],
                    &mut err,
                );
                errmax[(nd - 1) as usize] = err;
            }

            // -------------- IF ERRMAX > EPSAPR, stop ---------------------
            if errmax[(nd - 1) as usize] > epsapr[(nd - 1) as usize] {
                *iercod = -1;
                nu = *ndguli;
                nv = *ndgvli;

                // --------- Otherwise, try to remove again the coeff -------
            } else {
                minu = 1.max(((*iordru << 1) + 1)).max(*ndminu);
                maxu = *ndguli;
                minv = 1.max(((*iordrv << 1) + 1)).max(*ndminv);
                maxv = *ndgvli;
                if maxu >= (*iordru + 1) << 1 && maxv >= (*iordrv + 1) << 1 {
                    let mut ndses_ = ndses;
                    let mut minu_ = minu;
                    let mut maxu_ = maxu;
                    let mut minv_ = minv;
                    let mut maxv_ = maxv;
                    let mut epsapr_nd = epsapr[(nd - 1) as usize];
                    let mut err_nd = errmax[(nd - 1) as usize];
                    mma2er2_(
                        ndjacu,
                        ndjacv,
                        &mut ndses_,
                        &mut minu_,
                        &mut maxu_,
                        &mut minv_,
                        &mut maxv_,
                        iordru,
                        iordrv,
                        xmaxju,
                        xmaxjv,
                        &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                        &mut epsapr_nd,
                        &mut vecerr[..],
                        &mut err_nd,
                        &mut nu,
                        &mut nv,
                    );
                    errmax[(nd - 1) as usize] = err_nd;
                } else {
                    nu = maxu;
                    nv = maxv;
                }
            }

            // --------------------- Calculate the average error ------------
            let nu1 = nu + 1;
            let nv1 = nv + 1;
            let mut ndses_ = ndses;
            let mut nu1_ = nu1;
            let mut nv1_ = nv1;
            let mut ndjacu_a = *ndjacu;
            let mut ndjacu_b = *ndjacu;
            let mut ndjacv_a = *ndjacv;
            let mut ndjacv_b = *ndjacv;
            let mut errmoy_nd = errmoy[(nd - 1) as usize];
            mma2moy_(
                &mut ndjacu_a,
                &mut ndjacv_a,
                &mut ndses_,
                &mut nu1_,
                &mut ndjacu_b,
                &mut nv1_,
                &mut ndjacv_b,
                iordru,
                iordrv,
                &patjac[((idim - 1) * patjac_dim2 * patjac_dim1) as usize..],
                &mut errmoy_nd,
            );
            errmoy[(nd - 1) as usize] = errmoy_nd;

            // --------------------- Set to 0.D0 the rejected coeffs --------
            for ii in idim..=idim + ndses - 1 {
                for jj in nv1..=*ndjacv {
                    for kk in nu1..=*ndjacu {
                        patjac[((kk - 1)
                            + ((jj - 1) + (ii - 1) * patjac_dim2) * patjac_dim1)
                            as usize] = 0.;
                    }
                }
            }

            // --------------- Return the nb of coeff of approximation ------
            *ndegpu = (*ndegpu).max(nu);
            *ndegpv = (*ndegpv).max(nv);
            idim += ndses;
        }
    }

    // ------------------------------ The end -------------------------------
    // OCCT L9999:
    sys_base::maermsg_("MMA2CE2", iercod);
    if ldbg {
        sys_base::mgsomsg_("MMA2CE2");
    }
}

/// OCCT mma2ce2_ cut-management labels L9100 / L9200 / L9300 (each may fall
/// to L9001 when NUMDEC is out of range).
fn itydec_from_numdec(
    numdec: &mut i32,
    itydec: &mut i32,
    iercod: &mut i32,
    label: i32,
    ldbg: bool,
) {
    // Each OCCT label first re-checks NUMDEC (goto L9001 when invalid).
    if *numdec <= 0 || *numdec > 5 {
        // OCCT L9001:
        *iercod = 1;
        sys_base::maermsg_("MMA2CE2", iercod);
        if ldbg {
            sys_base::mgsomsg_("MMA2CE2");
        }
        return;
    }
    match label {
        // L9100: cut by U if possible.
        9100 => {
            if *numdec != 2 {
                *itydec = 1;
            } else {
                *itydec = 2;
            }
        }
        // L9200: cut by V if possible.
        9200 => {
            if *numdec != 1 {
                *itydec = 2;
            } else {
                *itydec = 1;
            }
        }
        // L9300: cut by both if possible.
        _ => {
            if *numdec == 5 {
                *itydec = 3;
            } else if *numdec == 2 || *numdec == 4 {
                *itydec = 2;
            } else if *numdec == 1 || *numdec == 3 {
                *itydec = 1;
            } else {
                // OCCT goto L9001.
                *iercod = 1;
                sys_base::maermsg_("MMA2CE2", iercod);
                if ldbg {
                    sys_base::mgsomsg_("MMA2CE2");
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// mma2ce1_ (AdvApp2Var_ApproxF2var.cxx L3719-4021)
// ---------------------------------------------------------------------------

/// OCCT mma2ce1_ (AdvApp2Var_ApproxF2var.cxx L3719-4021) - calculation of
/// coefficients of polynomial approximation of degree (NDJACU,NDJACV),
/// driver allocating the Gauss / Jacobi-maximum workspaces.  All buffers
/// are raw physical slices.
#[allow(clippy::too_many_arguments)]
pub fn mma2ce1_(
    numdec: &mut i32,
    ndimen: &mut i32,
    nbsesp: &mut i32,
    ndimse: &[i32],
    ndminu: &mut i32,
    ndminv: &mut i32,
    ndguli: &mut i32,
    ndgvli: &mut i32,
    ndjacu: &mut i32,
    ndjacv: &mut i32,
    iordru: &mut i32,
    iordrv: &mut i32,
    nbpntu: &mut i32,
    nbpntv: &mut i32,
    epsapr: &[f64],
    sosotb: &mut [f64],
    disotb: &mut [f64],
    soditb: &mut [f64],
    diditb: &mut [f64],
    patjac: &mut [f64],
    errmax: &mut [f64],
    errmoy: &mut [f64],
    ndegpu: &mut i32,
    ndegpv: &mut i32,
    itydec: &mut i32,
    iercod: &mut i32,
) {
    let c__8: i32 = 8;

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2CE1");
    }
    *iercod = 0;

    let isz1 = (*nbpntu / 2 + 1) * (*ndjacu - ((*iordru + 1) << 1) + 1);
    let isz2 = (*nbpntv / 2 + 1) * (*ndjacv - ((*iordrv + 1) << 1) + 1);
    let isz3 = (*nbpntv / 2 + 1) * (*ndjacu - ((*iordru + 1) << 1) + 1) * *ndimen;
    let isz4 = *nbpntv / 2 * (*ndjacu - ((*iordru + 1) << 1) + 1) * *ndimen;
    let isz5 = *ndjacu + 1 - ((*iordru + 1) << 1);
    let isz6 = *ndjacv + 1 - ((*iordrv + 1) << 1);
    let isz7 = *ndimen << 2;
    let mut iszwr = isz1 + isz2 + isz3 + isz4 + isz5 + isz6 + isz7;
    let mut wrkar: Vec<f64> = Vec::new();
    let mut iofwr: isize = 0;
    let mut ier = 0;
    let mut iunit = c__8;
    sys_base::mcrrqst_(&mut iunit, &mut iszwr, &mut wrkar, &mut iofwr, &mut ier);
    if ier > 0 {
        // OCCT goto L9013.
        *iercod = 13;
        sys_base::maermsg_("MMA2CE1", iercod);
        if ldbg {
            sys_base::mgsomsg_("MMA2CE1");
        }
        return;
    }
    let ipt1 = isz1;
    let ipt2 = ipt1 + isz2;
    let ipt3 = ipt2 + isz3;
    let ipt4 = ipt3 + isz4;
    let ipt5 = ipt4 + isz5;
    let ipt6 = ipt5 + isz6;

    {
        // Split the workspace into the OCCT zones:
        // [0..ipt1)=gssutb, [ipt1..ipt2)=gssvtb, [ipt2..ipt3)=chpair,
        // [ipt3..ipt4)=chimpr, [ipt4..ipt5)=xmaxju, [ipt5..ipt6)=xmaxjv,
        // [ipt6..iszwr)=vecerr.
        let (w_gssutb, w_rest) = wrkar.split_at_mut(ipt1 as usize);
        let (w_gssvtb, w_rest) = w_rest.split_at_mut(ipt2 as usize - ipt1 as usize);
        let (w_chpair, w_rest) = w_rest.split_at_mut(ipt3 as usize - ipt2 as usize);
        let (w_chimpr, w_rest) = w_rest.split_at_mut(ipt4 as usize - ipt3 as usize);
        let (w_xmaxju, w_rest) = w_rest.split_at_mut(ipt5 as usize - ipt4 as usize);
        let (w_xmaxjv, w_vecerr) = w_rest.split_at_mut(ipt6 as usize - ipt5 as usize);
        let _ = w_vecerr;

        // ----------------- Return Gauss coefficients of integration -----
        // OCCT: mmapptt_(ndjacu, nbpntu, iordru, wrkar_off, iercod).
        mmapptt_(ndjacu, nbpntu, iordru, w_gssutb, iercod);
        if *iercod > 0 {
            finish_mma2ce1_(ldbg, iercod);
            return;
        }
        // OCCT: mmapptt_(ndjacv, nbpntv, iordrv, &wrkar_off[ipt1], iercod).
        mmapptt_(ndjacv, nbpntv, iordrv, w_gssvtb, iercod);
        if *iercod > 0 {
            finish_mma2ce1_(ldbg, iercod);
            return;
        }

        // ------------------- Return max polynoms of  Jacobi --------------
        // OCCT (L3952-3953): both results are written to &wrkar_off[ipt5]
        // (the XMAXJV zone); the XMAXJU zone is never initialized here.
        let mut ndjacu_ = *ndjacu;
        let mut iordru_ = *iordru;
        mma2jmx_(&mut ndjacu_, &mut iordru_, w_xmaxjv);
        let mut ndjacv_ = *ndjacv;
        let mut iordrv_ = *iordrv;
        mma2jmx_(&mut ndjacv_, &mut iordrv_, w_xmaxjv);

        // ------ Calculate the coefficients and their error contribution --
        mma2ce2_(
            numdec,
            ndimen,
            nbsesp,
            ndimse,
            ndminu,
            ndminv,
            ndguli,
            ndgvli,
            ndjacu,
            ndjacv,
            iordru,
            iordrv,
            nbpntu,
            nbpntv,
            epsapr,
            sosotb,
            disotb,
            soditb,
            diditb,
            w_gssutb,
            w_gssvtb,
            w_xmaxju,
            w_xmaxjv,
            w_vecerr,
            w_chpair,
            w_chimpr,
            patjac,
            errmax,
            errmoy,
            ndegpu,
            ndegpv,
            itydec,
            iercod,
        );
    }
    if *iercod > 0 {
        finish_mma2ce1_(ldbg, iercod);
        return;
    }
    finish_mma2ce1_(ldbg, iercod);
}

/// OCCT mma2ce1_ L9999 tail (mcrdelt_ + maermsg_ + mgsomsg_).
fn finish_mma2ce1_(ldbg: bool, iercod: &mut i32) {
    let c__8: i32 = 8;
    // OCCT: if (iofwr != 0) mcrdelt_(...) - the rcad encoding always
    // releases the block after a successful allocation.
    {
        let mut ier = 0;
        let mut isize_ = 0;
        let mut iofst: isize = 0;
        let mut t: Vec<f64> = Vec::new();
        let mut iunit = c__8;
        sys_base::mcrdelt_(&mut iunit, &mut isize_, &mut t, &mut iofst, &mut ier);
        if ier > 0 {
            *iercod = 13;
        }
    }
    sys_base::maermsg_("MMA2CE1", iercod);
    if ldbg {
        sys_base::mgsomsg_("MMA2CE1");
    }
}
