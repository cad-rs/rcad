//! OCCT AdvApp2Var_ApproxF2var (AdvApp2Var_ApproxF2var.cxx) - part D: the
//! constraint-discretization engines (mma2can_, mma2cd1_, mma2cd2_,
//! mma2cd3_, mma2cdi_) consumed by AdvApp2Var_Patch::AddConstraints.
//!
//! Encoding note: as in parts A-C, the OCCT Fortran 1-based index formulas
//! are kept on 0-based slices by subtracting the OCCT virtual origin inside
//! the flat-index expression; `&arr[k]` (1-based k) passes as
//! `&arr[(k - 1) as usize..]` (or the raw sub-slice when the base was
//! already offset-adjusted).

use super::approxf2var_a::mma1her_;
use super::approxf2var_c::mmjacpt_;
use super::math_base_b::mmfmca8_;
use super::math_base_b::mmmpocur_;
use super::sys_base;

// ---------------------------------------------------------------------------
// mma2can_ (AdvApp2Var_ApproxF2var.cxx L2205-2363)
// ---------------------------------------------------------------------------

/// OCCT mma2can_ (AdvApp2Var_ApproxF2var.cxx L2205-2363) - change of Jacobi
/// base to canonical (-1,1) and writing in a greater table.  `patjac` /
/// `pataux` / `patcan` are raw physical slices (pataux holds
/// 2 * ncoefu * ncoefv * ndimen doubles).
pub fn mma2can_(
    ncfmxu: &i32,
    ncfmxv: &i32,
    ndimen: &i32,
    iordru: &i32,
    iordrv: &i32,
    ncoefu: &i32,
    ncoefv: &i32,
    patjac: &[f64],
    pataux: &mut [f64],
    patcan: &mut [f64],
    iercod: &mut i32,
) {
    // OCCT: patcan_dim1 = *ncfmxu; patcan_dim2 = *ncfmxv;
    //       patcan_offset = patcan_dim1 * (patcan_dim2 + 1) + 1.
    let patcan_dim1 = *ncfmxu;
    let patcan_dim2 = *ncfmxv;

    let ldbg = sys_base::mnfndeb_() >= 2;
    if ldbg {
        sys_base::mgenmsg_("MMA2CAN");
    }
    *iercod = 0;

    if *iordru < -1 || *iordru > 2 {
        // OCCT goto L9100.
        *iercod = 1;
        sys_base::maermsg_("MMA2CAN", iercod);
        if ldbg {
            sys_base::mgsomsg_("MMA2CAN");
        }
        return;
    }
    if *iordrv < -1 || *iordrv > 2 {
        *iercod = 1;
        sys_base::maermsg_("MMA2CAN", iercod);
        if ldbg {
            sys_base::mgsomsg_("MMA2CAN");
        }
        return;
    }
    if *ncoefu > *ncfmxu || *ncoefv > *ncfmxv {
        *iercod = 1;
        sys_base::maermsg_("MMA2CAN", iercod);
        if ldbg {
            sys_base::mgsomsg_("MMA2CAN");
        }
        return;
    }

    // --> Pass to canonic base (-1,1)
    // OCCT: mmjacpt_(ndimen, ncoefu, ncoefv, iordru, iordrv,
    //   &patjac[patjac_offset], &pataux[1], &patcan[patcan_offset])
    //   (all sub-pointers resolve to the raw addresses).
    mmjacpt_(ndimen, ncoefu, ncoefv, iordru, iordrv, patjac, pataux, patcan);

    // --> Write all in a greater table
    // OCCT: mmfmca8_(ncoefu, ncoefv, ndimen, ncfmxu, ncfmxv, ndimen,
    //   &patcan[patcan_offset], &patcan[patcan_offset]) - IN-PLACE call.
    // The expansion is safe in OCCT because the copies run backwards over
    // strictly increasing destinations; the rcad encoding snapshots the
    // read-only source (tabini is never written by mmfmca8_), which is
    // observationally equivalent.
    let tabini_snapshot: Vec<f64> = patcan.to_vec();
    mmfmca8_(
        ncoefu,
        ncoefv,
        ndimen,
        ncfmxu,
        ncfmxv,
        ndimen,
        &tabini_snapshot,
        patcan,
    );

    // --> Complete with zeros the resulting table.
    let ilon1 = *ncfmxu - *ncoefu;
    let ilon2 = *ncfmxu * (*ncfmxv - *ncoefv);
    for nd in 1..=*ndimen {
        if ilon1 > 0 {
            for ii in 1..=*ncoefv {
                // OCCT: mvriraz_(&ilon1,
                //   &patcan[*ncoefu + 1 + (ii + nd * dim2) * dim1])
                //   raw = (ncoefu - 1) + (ii - 1) * dim1 + (nd - 1) * dim1 * dim2.
                let raw = ((ncoefu - 1)
                    + (ii - 1) * patcan_dim1
                    + (nd - 1) * patcan_dim1 * patcan_dim2) as usize;
                sys_base::mvriraz_(ilon1, &mut patcan[raw..]);
            }
        }
        if ilon2 > 0 {
            // OCCT: mvriraz_(&ilon2,
            //   &patcan[(*ncoefv + 1 + nd * dim2) * dim1 + 1])
            //   raw = (ncoefv + (nd - 1) * dim2) * dim1.
            let raw = ((ncoefv + (nd - 1) * patcan_dim2) * patcan_dim1) as usize;
            sys_base::mvriraz_(ilon2, &mut patcan[raw..]);
        }
    }

    // ------------------------------ The end -------------------------------
    // OCCT L9999:
    sys_base::maermsg_("MMA2CAN", iercod);
    if ldbg {
        sys_base::mgsomsg_("MMA2CAN");
    }
}

// ---------------------------------------------------------------------------
// mma2cd1_ (AdvApp2Var_ApproxF2var.cxx L2365-2722)
// ---------------------------------------------------------------------------

/// OCCT mma2cd1_ (AdvApp2Var_ApproxF2var.cxx L2365-2722) - discretization of
/// the corner constraint interpolation polynoms of order IORDRE.  `urootl` /
/// `vrootl` are 1-based emulated; the other buffers are raw physical slices.
#[allow(clippy::too_many_arguments)]
pub fn mma2cd1_(
    ndimen: &i32,
    nbpntu: &i32,
    urootl: &[f64],
    nbpntv: &i32,
    vrootl: &[f64],
    iordru: &i32,
    iordrv: &i32,
    contr1: &[f64],
    contr2: &[f64],
    contr3: &[f64],
    contr4: &[f64],
    fpntbu: &mut [f64],
    fpntbv: &mut [f64],
    uhermt: &mut [f64],
    vhermt: &mut [f64],
    sosotb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
    diditb: &mut [f64],
) {
    let mut c__1: i32 = 1;

    // OCCT: sosotb_dim1 = *nbpntu / 2 + 1; sosotb_dim2 = *nbpntv / 2 + 1.
    let sosotb_dim1 = (nbpntu / 2 + 1) as i32;
    let sosotb_dim2 = (nbpntv / 2 + 1) as i32;
    // OCCT: diditb_dim1 = *nbpntu / 2 + 1; diditb_dim2 = *nbpntv / 2 + 1.
    let diditb_dim1 = (nbpntu / 2 + 1) as i32;
    let diditb_dim2 = (nbpntv / 2 + 1) as i32;
    // OCCT: soditb_dim1 = *nbpntu / 2; soditb_dim2 = *nbpntv / 2.
    let soditb_dim1 = (nbpntu / 2) as i32;
    let soditb_dim2 = (nbpntv / 2) as i32;
    // OCCT: uhermt_dim1 = (*iordru << 1) + 2; fpntbu_dim1 = *nbpntu.
    let uhermt_dim1 = ((*iordru << 1) + 2) as i32;
    let fpntbu_dim1 = *nbpntu;
    // OCCT: vhermt_dim1 = (*iordrv << 1) + 2; fpntbv_dim1 = *nbpntv.
    let vhermt_dim1 = ((*iordrv << 1) + 2) as i32;
    let fpntbv_dim1 = *nbpntv;
    // OCCT: contrX_dim1 = *ndimen; contrX_dim2 = *iordru + 2.
    let contr1_dim1 = *ndimen;
    let contr1_dim2 = *iordru + 2;

    let fpbu = |ll: i32, col: i32| -> usize {
        // OCCT: fpntbu[ll + col * fpntbu_dim1] - (dim1 + 1).
        ((ll - 1) + (col - 1) * fpntbu_dim1) as usize
    };
    let fpbv = |kk: i32, col: i32| -> usize {
        ((kk - 1) + (col - 1) * fpntbv_dim1) as usize
    };
    let ss = |ll: i32, kk: i32, nd: i32| -> usize {
        // OCCT: sosotb[ll + (kk + nd * dim2) * dim1] - dim1 * dim2.
        ((ll - 1) + ((kk - 1) + (nd - 1) * sosotb_dim2) * sosotb_dim1) as usize
    };
    let sd = |ll: i32, kk: i32, nd: i32| -> usize {
        // OCCT: soditb[ll + (kk + nd * dim2) * dim1] - dim1 * (dim2 + 1) - 1.
        ((ll - 1) + ((kk - 1) + (nd - 1) * soditb_dim2) * soditb_dim1) as usize
    };
    let dd = |ll: i32, kk: i32, nd: i32| -> usize {
        // OCCT: diditb[ll + (kk + nd * dim2) * dim1] - dim1 * dim2.
        ((ll - 1) + ((kk - 1) + (nd - 1) * diditb_dim2) * diditb_dim1) as usize
    };
    let ct = |nd: i32, ii: i32, jj: i32| -> usize {
        // OCCT: contrX[nd + (ii + jj * dim2) * dim1] - dim1 * (dim2 + 1) - 1.
        ((nd - 1) + ((ii - 1) + (jj - 1) * contr1_dim2) * contr1_dim1) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA2CD1");
    }

    // ------------------- Discretisation of Hermite polynoms -----------
    let ncfhu = (*iordru + 1) << 1;
    for ii in 1..=ncfhu {
        for ll in 1..=*nbpntu {
            // OCCT: mmmpocur_(&ncfhu, &c__1, &ncfhu, &uhermt[ii * uhermt_dim1],
            //   &urootl[ll], &fpntbu[ll + ii * fpntbu_dim1]).
            // OCCT aliases the same int* for NCOFMX/NDEG; the callee only
            // reads them, so the rcad encoding splits the alias.
            let mut ncfhu_a = ncfhu;
            let mut ncfhu_b = ncfhu;
            let mut tparam = urootl[(ll - 1) as usize];
            mmmpocur_(
                &mut ncfhu_a,
                &mut c__1,
                &mut ncfhu_b,
                &uhermt[((ii - 1) * uhermt_dim1) as usize..],
                &mut tparam,
                &mut fpntbu[fpbu(ll, ii)..],
            );
        }
    }
    let ncfhv = (*iordrv + 1) << 1;
    for jj in 1..=ncfhv {
        for kk in 1..=*nbpntv {
            let mut ncfhv_a = ncfhv;
            let mut ncfhv_b = ncfhv;
            let mut tparam = vrootl[(kk - 1) as usize];
            mmmpocur_(
                &mut ncfhv_a,
                &mut c__1,
                &mut ncfhv_b,
                &vhermt[((jj - 1) * vhermt_dim1) as usize..],
                &mut tparam,
                &mut fpntbv[fpbv(kk, jj)..],
            );
        }
    }

    // ---- The discretizations of polynoms of constraints are subtracted ----
    let nuroo = *nbpntu / 2;
    let nvroo = *nbpntv / 2;
    for nd in 1..=*ndimen {
        for jj in 1..=*iordrv + 1 {
            for ii in 1..=*iordru + 1 {
                let bid1 = contr1[ct(nd, ii, jj)];
                let bid2 = contr2[ct(nd, ii, jj)];
                let bid3 = contr3[ct(nd, ii, jj)];
                let bid4 = contr4[ct(nd, ii, jj)];

                for kk in 1..=nvroo {
                    let kkp = (*nbpntv + 1) / 2 + kk;
                    let kkm = nvroo - kk + 1;
                    let sov1 = fpntbv[fpbv(kkp, (jj << 1) - 1)] + fpntbv[fpbv(kkm, (jj << 1) - 1)];
                    let div1 = fpntbv[fpbv(kkp, (jj << 1) - 1)] - fpntbv[fpbv(kkm, (jj << 1) - 1)];
                    let sov2 = fpntbv[fpbv(kkp, jj << 1)] + fpntbv[fpbv(kkm, jj << 1)];
                    let div2 = fpntbv[fpbv(kkp, jj << 1)] - fpntbv[fpbv(kkm, jj << 1)];
                    for ll in 1..=nuroo {
                        let llp = (*nbpntu + 1) / 2 + ll;
                        let llm = nuroo - ll + 1;
                        let sou1 =
                            fpntbu[fpbu(llp, (ii << 1) - 1)] + fpntbu[fpbu(llm, (ii << 1) - 1)];
                        let diu1 =
                            fpntbu[fpbu(llp, (ii << 1) - 1)] - fpntbu[fpbu(llm, (ii << 1) - 1)];
                        let sou2 = fpntbu[fpbu(llp, ii << 1)] + fpntbu[fpbu(llm, ii << 1)];
                        let diu2 = fpntbu[fpbu(llp, ii << 1)] - fpntbu[fpbu(llm, ii << 1)];
                        sosotb[ss(ll, kk, nd)] = sosotb[ss(ll, kk, nd)]
                            - bid1 * sou1 * sov1
                            - bid2 * sou2 * sov1
                            - bid3 * sou1 * sov2
                            - bid4 * sou2 * sov2;
                        soditb[sd(ll, kk, nd)] = soditb[sd(ll, kk, nd)]
                            - bid1 * sou1 * div1
                            - bid2 * sou2 * div1
                            - bid3 * sou1 * div2
                            - bid4 * sou2 * div2;
                        disotb[sd(ll, kk, nd)] = disotb[sd(ll, kk, nd)]
                            - bid1 * diu1 * sov1
                            - bid2 * diu2 * sov1
                            - bid3 * diu1 * sov2
                            - bid4 * diu2 * sov2;
                        diditb[dd(ll, kk, nd)] = diditb[dd(ll, kk, nd)]
                            - bid1 * diu1 * div1
                            - bid2 * diu2 * div1
                            - bid3 * diu1 * div2
                            - bid4 * diu2 * div2;
                    }
                }

                // ------------ Case when the discretization is done only on
                // the roots of Legendre polynom of uneven degree, 0 is root
                if *nbpntu % 2 == 1 {
                    let sou1 = fpntbu[fpbu(nuroo + 1, (ii << 1) - 1)];
                    let sou2 = fpntbu[fpbu(nuroo + 1, ii << 1)];
                    for kk in 1..=nvroo {
                        let kkp = (*nbpntv + 1) / 2 + kk;
                        let kkm = nvroo - kk + 1;
                        let sov1 =
                            fpntbv[fpbv(kkp, (jj << 1) - 1)] + fpntbv[fpbv(kkm, (jj << 1) - 1)];
                        let div1 =
                            fpntbv[fpbv(kkp, (jj << 1) - 1)] - fpntbv[fpbv(kkm, (jj << 1) - 1)];
                        let sov2 = fpntbv[fpbv(kkp, jj << 1)] + fpntbv[fpbv(kkm, jj << 1)];
                        let div2 = fpntbv[fpbv(kkp, jj << 1)] - fpntbv[fpbv(kkm, jj << 1)];
                        // OCCT: sosotb[(kk + nd * dim2) * dim1] (ll slot 0).
                        sosotb[ss(0, kk, nd)] = sosotb[ss(0, kk, nd)]
                            - bid1 * sou1 * sov1
                            - bid2 * sou2 * sov1
                            - bid3 * sou1 * sov2
                            - bid4 * sou2 * sov2;
                        diditb[dd(0, kk, nd)] = diditb[dd(0, kk, nd)]
                            - bid1 * sou1 * div1
                            - bid2 * sou2 * div1
                            - bid3 * sou1 * div2
                            - bid4 * sou2 * div2;
                    }
                }

                if *nbpntv % 2 == 1 {
                    let sov1 = fpntbv[fpbv(nvroo + 1, (jj << 1) - 1)];
                    let sov2 = fpntbv[fpbv(nvroo + 1, jj << 1)];
                    for ll in 1..=nuroo {
                        let llp = (*nbpntu + 1) / 2 + ll;
                        let llm = nuroo - ll + 1;
                        let sou1 =
                            fpntbu[fpbu(llp, (ii << 1) - 1)] + fpntbu[fpbu(llm, (ii << 1) - 1)];
                        let diu1 =
                            fpntbu[fpbu(llp, (ii << 1) - 1)] - fpntbu[fpbu(llm, (ii << 1) - 1)];
                        let sou2 = fpntbu[fpbu(llp, ii << 1)] + fpntbu[fpbu(llm, ii << 1)];
                        let diu2 = fpntbu[fpbu(llp, ii << 1)] - fpntbu[fpbu(llm, ii << 1)];
                        // OCCT: sosotb[ll + nd * dim2 * dim1] (kk slot 0).
                        sosotb[ss(ll, 0, nd)] = sosotb[ss(ll, 0, nd)]
                            - bid1 * sou1 * sov1
                            - bid2 * sou2 * sov1
                            - bid3 * sou1 * sov2
                            - bid4 * sou2 * sov2;
                        diditb[dd(ll, 0, nd)] = diditb[dd(ll, 0, nd)]
                            - bid1 * diu1 * sov1
                            - bid2 * diu2 * sov1
                            - bid3 * diu1 * sov2
                            - bid4 * diu2 * sov2;
                    }
                }

                if *nbpntu % 2 == 1 && *nbpntv % 2 == 1 {
                    let sou1 = fpntbu[fpbu(nuroo + 1, (ii << 1) - 1)];
                    let sou2 = fpntbu[fpbu(nuroo + 1, ii << 1)];
                    let sov1 = fpntbv[fpbv(nvroo + 1, (jj << 1) - 1)];
                    let sov2 = fpntbv[fpbv(nvroo + 1, jj << 1)];
                    sosotb[ss(0, 0, nd)] = sosotb[ss(0, 0, nd)]
                        - bid1 * sou1 * sov1
                        - bid2 * sou2 * sov1
                        - bid3 * sou1 * sov2
                        - bid4 * sou2 * sov2;
                    diditb[dd(0, 0, nd)] = diditb[dd(0, 0, nd)]
                        - bid1 * sou1 * sov1
                        - bid2 * sou2 * sov1
                        - bid3 * sou1 * sov2
                        - bid4 * sou2 * sov2;
                }
            }
        }
    }

    // ------------------------------ The End -------------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA2CD1");
    }
}

// ---------------------------------------------------------------------------
// mma2cd2_ (AdvApp2Var_ApproxF2var.cxx L2723-3041)
// ---------------------------------------------------------------------------

/// OCCT mma2cd2_ (AdvApp2Var_ApproxF2var.cxx L2723-3041) - discretization of
/// the border constraint polynoms on the 2 iso-V borders of order IORDRV.
#[allow(clippy::too_many_arguments)]
pub fn mma2cd2_(
    ndimen: &i32,
    nbpntu: &i32,
    nbpntv: &i32,
    vrootl: &[f64],
    iordrv: &i32,
    sotbv1: &[f64],
    sotbv2: &[f64],
    ditbv1: &[f64],
    ditbv2: &[f64],
    fpntab: &mut [f64],
    vhermt: &mut [f64],
    sosotb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
    diditb: &mut [f64],
) {
    let mut c__1: i32 = 1;

    let sosotb_dim1 = (nbpntu / 2 + 1) as i32;
    let sosotb_dim2 = (nbpntv / 2 + 1) as i32;
    let diditb_dim1 = (nbpntu / 2 + 1) as i32;
    let diditb_dim2 = (nbpntv / 2 + 1) as i32;
    let soditb_dim1 = (nbpntu / 2) as i32;
    let soditb_dim2 = (nbpntv / 2) as i32;
    let vhermt_dim1 = ((*iordrv << 1) + 2) as i32;
    // OCCT: fpntab_dim1 = *nbpntv.
    let fpntab_dim1 = *nbpntv;
    // OCCT: sotbvX_dim1 = *nbpntu / 2 + 1; sotbvX_dim2 = *ndimen.
    let sotbv_dim1 = (nbpntu / 2 + 1) as i32;
    let sotbv_dim2 = *ndimen;

    let fpa = |jj: i32, col: i32| -> usize {
        // OCCT: fpntab[jj + col * fpntab_dim1] - (dim1 + 1).
        ((jj - 1) + (col - 1) * fpntab_dim1) as usize
    };
    let ss = |kk: i32, jj: i32, nd: i32| -> usize {
        // OCCT: sosotb[kk + (jj + nd * dim2) * dim1] - dim1 * dim2.
        ((kk - 1) + ((jj - 1) + (nd - 1) * sosotb_dim2) * sosotb_dim1) as usize
    };
    let sd = |kk: i32, jj: i32, nd: i32| -> usize {
        ((kk - 1) + ((jj - 1) + (nd - 1) * soditb_dim2) * soditb_dim1) as usize
    };
    let dd = |kk: i32, jj: i32, nd: i32| -> usize {
        ((kk - 1) + ((jj - 1) + (nd - 1) * diditb_dim2) * diditb_dim1) as usize
    };
    let sv = |kk: i32, nd: i32, ii: i32| -> usize {
        // OCCT: sotbvX[kk + (nd + ii * dim2) * dim1] - dim1 * (dim2 + 1).
        ((kk - 1) + ((nd - 1) + (ii - 1) * sotbv_dim2) * sotbv_dim1) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA2CD2");
    }

    // ------------------- Discretization of Hermit polynoms -----------
    let ncfhv = (*iordrv + 1) << 1;
    for ii in 1..=ncfhv {
        for jj in 1..=*nbpntv {
            let mut ncfhv_a = ncfhv;
            let mut ncfhv_b = ncfhv;
            let mut tparam = vrootl[(jj - 1) as usize];
            mmmpocur_(
                &mut ncfhv_a,
                &mut c__1,
                &mut ncfhv_b,
                &vhermt[((ii - 1) * vhermt_dim1) as usize..],
                &mut tparam,
                &mut fpntab[fpa(jj, ii)..],
            );
        }
    }

    // ---- The discretizations of polynoms of constraints are subtracted ----
    let nuroo = *nbpntu / 2;
    let nvroo = *nbpntv / 2;

    for nd in 1..=*ndimen {
        for ii in 1..=*iordrv + 1 {
            for kk in 1..=nuroo {
                let bid1 = sotbv1[sv(kk, nd, ii)];
                let bid2 = sotbv2[sv(kk, nd, ii)];
                let bid3 = ditbv1[sv(kk, nd, ii)];
                let bid4 = ditbv2[sv(kk, nd, ii)];
                for jj in 1..=nvroo {
                    let jjp = (*nbpntv + 1) / 2 + jj;
                    let jjm = nvroo - jj + 1;
                    sosotb[ss(kk, jj, nd)] = sosotb[ss(kk, jj, nd)]
                        - bid1
                            * (fpntab[fpa(jjp, (ii << 1) - 1)]
                                + fpntab[fpa(jjm, (ii << 1) - 1)])
                        - bid2 * (fpntab[fpa(jjp, ii << 1)] + fpntab[fpa(jjm, ii << 1)]);
                    disotb[sd(kk, jj, nd)] = disotb[sd(kk, jj, nd)]
                        - bid3
                            * (fpntab[fpa(jjp, (ii << 1) - 1)]
                                + fpntab[fpa(jjm, (ii << 1) - 1)])
                        - bid4 * (fpntab[fpa(jjp, ii << 1)] + fpntab[fpa(jjm, ii << 1)]);
                    soditb[sd(kk, jj, nd)] = soditb[sd(kk, jj, nd)]
                        - bid1
                            * (fpntab[fpa(jjp, (ii << 1) - 1)]
                                - fpntab[fpa(jjm, (ii << 1) - 1)])
                        - bid2 * (fpntab[fpa(jjp, ii << 1)] - fpntab[fpa(jjm, ii << 1)]);
                    diditb[dd(kk, jj, nd)] = diditb[dd(kk, jj, nd)]
                        - bid3
                            * (fpntab[fpa(jjp, (ii << 1) - 1)]
                                - fpntab[fpa(jjm, (ii << 1) - 1)])
                        - bid4 * (fpntab[fpa(jjp, ii << 1)] - fpntab[fpa(jjm, ii << 1)]);
                }
            }
        }

        // ------------ Case when the discretization is done only on the roots
        // of Legendre polynom of uneven degree, 0 is root
        if *nbpntv % 2 == 1 {
            for ii in 1..=*iordrv + 1 {
                for kk in 1..=nuroo {
                    // OCCT: sosotb[kk + nd * dim2 * dim1] (jj slot 0).
                    let bid1 = sotbv1[sv(kk, nd, ii)]
                        * fpntab[fpa(nvroo + 1, (ii << 1) - 1)]
                        + sotbv2[sv(kk, nd, ii)] * fpntab[fpa(nvroo + 1, ii << 1)];
                    sosotb[ss(kk, 0, nd)] -= bid1;
                    let bid2 = ditbv1[sv(kk, nd, ii)]
                        * fpntab[fpa(nvroo + 1, (ii << 1) - 1)]
                        + ditbv2[sv(kk, nd, ii)] * fpntab[fpa(nvroo + 1, ii << 1)];
                    diditb[dd(kk, 0, nd)] -= bid2;
                }
            }
        }

        if *nbpntu % 2 == 1 {
            for ii in 1..=*iordrv + 1 {
                for jj in 1..=nvroo {
                    let jjp = (*nbpntv + 1) / 2 + jj;
                    let jjm = nvroo - jj + 1;
                    // OCCT: sotbvX[(nd + ii * dim2) * dim1] (kk slot 0);
                    // sosotb[(jj + nd * dim2) * dim1] (kk slot 0).
                    let bid1 = sotbv1[sv(0, nd, ii)]
                        * (fpntab[fpa(jjp, (ii << 1) - 1)] + fpntab[fpa(jjm, (ii << 1) - 1)])
                        + sotbv2[sv(0, nd, ii)]
                            * (fpntab[fpa(jjp, ii << 1)] + fpntab[fpa(jjm, ii << 1)]);
                    sosotb[ss(0, jj, nd)] -= bid1;
                    let bid2 = sotbv1[sv(0, nd, ii)]
                        * (fpntab[fpa(jjp, (ii << 1) - 1)] - fpntab[fpa(jjm, (ii << 1) - 1)])
                        + sotbv2[sv(0, nd, ii)]
                            * (fpntab[fpa(jjp, ii << 1)] - fpntab[fpa(jjm, ii << 1)]);
                    diditb[dd(0, jj, nd)] -= bid2;
                }
            }
        }

        if *nbpntu % 2 == 1 && *nbpntv % 2 == 1 {
            for ii in 1..=*iordrv + 1 {
                let bid1 = sotbv1[sv(0, nd, ii)]
                    * fpntab[fpa(nvroo + 1, (ii << 1) - 1)]
                    + sotbv2[sv(0, nd, ii)] * fpntab[fpa(nvroo + 1, ii << 1)];
                sosotb[ss(0, 0, nd)] -= bid1;
            }
        }
    }

    // ------------------------------ The End -------------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA2CD2");
    }
}

// ---------------------------------------------------------------------------
// mma2cd3_ (AdvApp2Var_ApproxF2var.cxx L3043-3364)
// ---------------------------------------------------------------------------

/// OCCT mma2cd3_ (AdvApp2Var_ApproxF2var.cxx L3043-3364) - discretization of
/// the border constraint polynoms on the 2 iso-U borders of order IORDRU.
#[allow(clippy::too_many_arguments)]
pub fn mma2cd3_(
    ndimen: &i32,
    nbpntu: &i32,
    urootl: &[f64],
    nbpntv: &i32,
    iordru: &i32,
    sotbu1: &[f64],
    sotbu2: &[f64],
    ditbu1: &[f64],
    ditbu2: &[f64],
    fpntab: &mut [f64],
    uhermt: &mut [f64],
    sosotb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
    diditb: &mut [f64],
) {
    let mut c__1: i32 = 1;

    let sosotb_dim1 = (nbpntu / 2 + 1) as i32;
    let sosotb_dim2 = (nbpntv / 2 + 1) as i32;
    let diditb_dim1 = (nbpntu / 2 + 1) as i32;
    let diditb_dim2 = (nbpntv / 2 + 1) as i32;
    let soditb_dim1 = (nbpntu / 2) as i32;
    let soditb_dim2 = (nbpntv / 2) as i32;
    let uhermt_dim1 = ((*iordru << 1) + 2) as i32;
    // OCCT: fpntab_dim1 = *nbpntu.
    let fpntab_dim1 = *nbpntu;
    // OCCT: sotbuX_dim1 = *nbpntv / 2 + 1; sotbuX_dim2 = *ndimen.
    let sotbu_dim1 = (nbpntv / 2 + 1) as i32;
    let sotbu_dim2 = *ndimen;

    let fpa = |kk: i32, col: i32| -> usize {
        ((kk - 1) + (col - 1) * fpntab_dim1) as usize
    };
    let ss = |kk: i32, jj: i32, nd: i32| -> usize {
        ((kk - 1) + ((jj - 1) + (nd - 1) * sosotb_dim2) * sosotb_dim1) as usize
    };
    let sd = |kk: i32, jj: i32, nd: i32| -> usize {
        ((kk - 1) + ((jj - 1) + (nd - 1) * soditb_dim2) * soditb_dim1) as usize
    };
    let dd = |kk: i32, jj: i32, nd: i32| -> usize {
        ((kk - 1) + ((jj - 1) + (nd - 1) * diditb_dim2) * diditb_dim1) as usize
    };
    let su = |jj: i32, nd: i32, ii: i32| -> usize {
        // OCCT: sotbuX[jj + (nd + ii * dim2) * dim1] - dim1 * (dim2 + 1).
        ((jj - 1) + ((nd - 1) + (ii - 1) * sotbu_dim2) * sotbu_dim1) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA2CD3");
    }

    // ------------------- Discretization of polynoms of Hermit -----------
    let ncfhu = (*iordru + 1) << 1;
    for ii in 1..=ncfhu {
        for kk in 1..=*nbpntu {
            let mut ncfhu_a = ncfhu;
            let mut ncfhu_b = ncfhu;
            let mut tparam = urootl[(kk - 1) as usize];
            mmmpocur_(
                &mut ncfhu_a,
                &mut c__1,
                &mut ncfhu_b,
                &uhermt[((ii - 1) * uhermt_dim1) as usize..],
                &mut tparam,
                &mut fpntab[fpa(kk, ii)..],
            );
        }
    }

    // ---- The discretizations of polynoms of constraints are subtracted ----
    let nvroo = *nbpntv / 2;
    let nuroo = *nbpntu / 2;

    for nd in 1..=*ndimen {
        for ii in 1..=*iordru + 1 {
            for jj in 1..=nvroo {
                let bid1 = sotbu1[su(jj, nd, ii)];
                let bid2 = sotbu2[su(jj, nd, ii)];
                let bid3 = ditbu1[su(jj, nd, ii)];
                let bid4 = ditbu2[su(jj, nd, ii)];
                for kk in 1..=nuroo {
                    let kkp = (*nbpntu + 1) / 2 + kk;
                    let kkm = nuroo - kk + 1;
                    sosotb[ss(kk, jj, nd)] = sosotb[ss(kk, jj, nd)]
                        - bid1
                            * (fpntab[fpa(kkp, (ii << 1) - 1)]
                                + fpntab[fpa(kkm, (ii << 1) - 1)])
                        - bid2 * (fpntab[fpa(kkp, ii << 1)] + fpntab[fpa(kkm, ii << 1)]);
                    disotb[sd(kk, jj, nd)] = disotb[sd(kk, jj, nd)]
                        - bid1
                            * (fpntab[fpa(kkp, (ii << 1) - 1)]
                                - fpntab[fpa(kkm, (ii << 1) - 1)])
                        - bid2 * (fpntab[fpa(kkp, ii << 1)] - fpntab[fpa(kkm, ii << 1)]);
                    soditb[sd(kk, jj, nd)] = soditb[sd(kk, jj, nd)]
                        - bid3
                            * (fpntab[fpa(kkp, (ii << 1) - 1)]
                                + fpntab[fpa(kkm, (ii << 1) - 1)])
                        - bid4 * (fpntab[fpa(kkp, ii << 1)] + fpntab[fpa(kkm, ii << 1)]);
                    diditb[dd(kk, jj, nd)] = diditb[dd(kk, jj, nd)]
                        - bid3
                            * (fpntab[fpa(kkp, (ii << 1) - 1)]
                                - fpntab[fpa(kkm, (ii << 1) - 1)])
                        - bid4 * (fpntab[fpa(kkp, ii << 1)] - fpntab[fpa(kkm, ii << 1)]);
                }
            }
        }

        // ------------ Case when the discretization is done only on the roots
        // of Legendre polynom of uneven degree, 0 is root
        if *nbpntu % 2 == 1 {
            for ii in 1..=*iordru + 1 {
                for jj in 1..=nvroo {
                    // OCCT: diditb[(jj + nd * dim2) * dim1] (kk slot 0).
                    let bid1 = sotbu1[su(jj, nd, ii)]
                        * fpntab[fpa(nuroo + 1, (ii << 1) - 1)]
                        + sotbu2[su(jj, nd, ii)] * fpntab[fpa(nuroo + 1, ii << 1)];
                    sosotb[ss(0, jj, nd)] -= bid1;
                    let bid2 = ditbu1[su(jj, nd, ii)]
                        * fpntab[fpa(nuroo + 1, (ii << 1) - 1)]
                        + ditbu2[su(jj, nd, ii)] * fpntab[fpa(nuroo + 1, ii << 1)];
                    diditb[dd(0, jj, nd)] -= bid2;
                }
            }
        }

        if *nbpntv % 2 == 1 {
            for ii in 1..=*iordru + 1 {
                for kk in 1..=nuroo {
                    let kkp = (*nbpntu + 1) / 2 + kk;
                    let kkm = nuroo - kk + 1;
                    // OCCT: sotbuX[(nd + ii * dim2) * dim1] (jj slot 0);
                    // sosotb[kk + nd * dim2 * dim1] (jj slot 0).
                    let bid1 = sotbu1[su(0, nd, ii)]
                        * (fpntab[fpa(kkp, (ii << 1) - 1)] + fpntab[fpa(kkm, (ii << 1) - 1)])
                        + sotbu2[su(0, nd, ii)]
                            * (fpntab[fpa(kkp, ii << 1)] + fpntab[fpa(kkm, ii << 1)]);
                    sosotb[ss(kk, 0, nd)] -= bid1;
                    let bid2 = sotbu1[su(0, nd, ii)]
                        * (fpntab[fpa(kkp, (ii << 1) - 1)] - fpntab[fpa(kkm, (ii << 1) - 1)])
                        + sotbu2[su(0, nd, ii)]
                            * (fpntab[fpa(kkp, ii << 1)] - fpntab[fpa(kkm, ii << 1)]);
                    diditb[dd(kk, 0, nd)] -= bid2;
                }
            }
        }

        if *nbpntu % 2 == 1 && *nbpntv % 2 == 1 {
            for ii in 1..=*iordru + 1 {
                let bid1 = sotbu1[su(0, nd, ii)]
                    * fpntab[fpa(nuroo + 1, (ii << 1) - 1)]
                    + sotbu2[su(0, nd, ii)] * fpntab[fpa(nuroo + 1, ii << 1)];
                sosotb[ss(0, 0, nd)] -= bid1;
            }
        }
    }

    // ------------------------------ The End -------------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA2CD3");
    }
}

// ---------------------------------------------------------------------------
// mma2cdi_ (AdvApp2Var_ApproxF2var.cxx L3366-3717)
// ---------------------------------------------------------------------------

/// OCCT mma2cdi_ (AdvApp2Var_ApproxF2var.cxx L3366-3717) - discretization of
/// the order-IORDRE constraint interpolation polynoms.  `urootl` / `vrootl`
/// are 1-based emulated; the constraint / sotb / sosotb-family buffers are
/// raw physical slices.
#[allow(clippy::too_many_arguments)]
pub fn mma2cdi_(
    ndimen: &mut i32,
    nbpntu: &mut i32,
    urootl: &[f64],
    nbpntv: &mut i32,
    vrootl: &[f64],
    iordru: &mut i32,
    iordrv: &mut i32,
    contr1: &[f64],
    contr2: &[f64],
    contr3: &[f64],
    contr4: &[f64],
    sotbu1: &[f64],
    sotbu2: &[f64],
    ditbu1: &[f64],
    ditbu2: &[f64],
    sotbv1: &[f64],
    sotbv2: &[f64],
    ditbv1: &[f64],
    ditbv2: &[f64],
    sosotb: &mut [f64],
    soditb: &mut [f64],
    disotb: &mut [f64],
    diditb: &mut [f64],
    iercod: &mut i32,
) {
    let c__8: i32 = 8;

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA2CDI");
    }
    *iercod = 0;
    if *iordru < -1 || *iordru > 2 {
        // OCCT goto L9100.
        *iercod = 1;
        sys_base::maermsg_("MMA2CDI", iercod);
        if ibb >= 3 {
            sys_base::mgsomsg_("MMA2CDI");
        }
        return;
    }
    if *iordrv < -1 || *iordrv > 2 {
        *iercod = 1;
        sys_base::maermsg_("MMA2CDI", iercod);
        if ibb >= 3 {
            sys_base::mgsomsg_("MMA2CDI");
        }
        return;
    }

    // ------------------------- Set to zero --------------------------------
    let mut ilong = (*nbpntu / 2 + 1) * (*nbpntv / 2 + 1) * *ndimen;
    sys_base::mvriraz_(ilong, sosotb);
    sys_base::mvriraz_(ilong, diditb);
    ilong = *nbpntu / 2 * (*nbpntv / 2) * *ndimen;
    sys_base::mvriraz_(ilong, soditb);
    sys_base::mvriraz_(ilong, disotb);
    if *iordru == -1 && *iordrv == -1 {
        // OCCT goto L9999 (nothing to do).
        sys_base::maermsg_("MMA2CDI", iercod);
        if ibb >= 3 {
            sys_base::mgsomsg_("MMA2CDI");
        }
        return;
    }

    let isz1 = ((*iordru + 1) << 2) * (*iordru + 1);
    let isz2 = ((*iordrv + 1) << 2) * (*iordrv + 1);
    let isz3 = ((*iordru + 1) << 1) * *nbpntu;
    let isz4 = ((*iordrv + 1) << 1) * *nbpntv;
    let mut iszwr = isz1 + isz2 + isz3 + isz4;
    let mut wrkar: Vec<f64> = Vec::new();
    let mut iofwr: isize = 0;
    let mut ier = 0;
    let mut iunit = c__8;
    sys_base::mcrrqst_(&mut iunit, &mut iszwr, &mut wrkar, &mut iofwr, &mut ier);
    if ier > 0 {
        // OCCT goto L9013.
        *iercod = 13;
        sys_base::maermsg_("MMA2CDI", iercod);
        if ibb >= 3 {
            sys_base::mgsomsg_("MMA2CDI");
        }
        return;
    }
    // OCCT: wrkar_off = the raw workspace base (iofwr collapses to 0 in the
    // rcad encoding); ipt1..ipt3 are 0-based sub-slice offsets.
    let ipt1 = isz1;
    let ipt2 = ipt1 + isz2;
    let ipt3 = ipt2 + isz3;

    if *iordru >= 0 && *iordru <= 2 {
        // --- Return 2*(IORDRU+1) coeff of 2*(IORDRU+1) polynoms of Hermite
        // OCCT: mma1her_(iordru, wrkar_off, iercod) (raw base).
        mma1her_(iordru, &mut wrkar[..], iercod);
        if *iercod > 0 {
            // OCCT goto L9100.
            *iercod = 1;
            sys_base::maermsg_("MMA2CDI", iercod);
            if ibb >= 3 {
                sys_base::mgsomsg_("MMA2CDI");
            }
            return;
        }

        // ---- Subtract discretizations of polynoms of constraints ----
        // OCCT: mma2cd3_(ndimen, nbpntu, &urootl[1], nbpntv, iordru,
        //   &sotbu1[1], &sotbu2[1], &ditbu1[1], &ditbu2[1],
        //   &wrkar_off[ipt2], wrkar_off, &sosotb[sosotb_offset],
        //   &soditb[soditb_offset], &disotb[disotb_offset],
        //   &diditb[diditb_offset]) (all sub-pointers raw addresses).
        // The rcad encoding splits the workspace into the disjoint
        // fpntbu (ipt2..ipt3) and uhermt (0..isz1) zones.
        {
            let (w_uhermt, w_rest) = wrkar.split_at_mut(ipt2 as usize);
            let (w_fpntbu, _) = w_rest.split_at_mut(ipt3 as usize - ipt2 as usize);
            mma2cd3_(
                ndimen,
                nbpntu,
                urootl,
                nbpntv,
                iordru,
                sotbu1,
                sotbu2,
                ditbu1,
                ditbu2,
                w_fpntbu,
                w_uhermt,
                sosotb,
                soditb,
                disotb,
                diditb,
            );
        }
    }

    if *iordrv >= 0 && *iordrv <= 2 {
        // --- Return 2*(IORDRV+1) coeff of 2*(IORDRV+1) polynoms of Hermite
        // OCCT: mma1her_(iordrv, &wrkar_off[ipt1], iercod).
        mma1her_(iordrv, &mut wrkar[ipt1 as usize..ipt2 as usize], iercod);
        if *iercod > 0 {
            *iercod = 1;
            sys_base::maermsg_("MMA2CDI", iercod);
            if ibb >= 3 {
                sys_base::mgsomsg_("MMA2CDI");
            }
            return;
        }

        // ---- Subtract discretisations of polynoms of constraint ----
        // OCCT: mma2cd2_(ndimen, nbpntu, nbpntv, &vrootl[1], iordrv,
        //   &sotbv1[1], &sotbv2[1], &ditbv1[1], &ditbv2[1],
        //   &wrkar_off[ipt3], &wrkar_off[ipt1], &sosotb[sosotb_offset],
        //   &soditb[soditb_offset], &disotb[disotb_offset],
        //   &diditb[diditb_offset]).
        // The rcad encoding splits the workspace into the disjoint
        // fpntab (ipt3..iszwr) and vhermt (ipt1..ipt2) zones.
        {
            let (_, w_rest) = wrkar.split_at_mut(ipt1 as usize);
            let (w_vhermt, w_rest2) = w_rest.split_at_mut(ipt2 as usize - ipt1 as usize);
            let (_, w_rest3) = w_rest2.split_at_mut(ipt3 as usize - ipt2 as usize);
            mma2cd2_(
                ndimen,
                nbpntu,
                nbpntv,
                vrootl,
                iordrv,
                sotbv1,
                sotbv2,
                ditbv1,
                ditbv2,
                w_rest3,
                w_vhermt,
                sosotb,
                soditb,
                disotb,
                diditb,
            );
        }
    }

    // --------------- Subtract constraints of corners ----------------
    if *iordru >= 0 && *iordrv >= 0 {
        // OCCT: mma2cd1_(ndimen, nbpntu, &urootl[1], nbpntv, &vrootl[1],
        //   iordru, iordrv, &contr1[contr1_offset], ..., &contr4[...],
        //   &wrkar_off[ipt2], &wrkar_off[ipt3], wrkar_off,
        //   &wrkar_off[ipt1], &sosotb[sosotb_offset], &soditb[...],
        //   &disotb[...], &diditb[diditb_offset]).
        // The rcad encoding splits the workspace into the disjoint
        // uhermt (0..isz1) / vhermt (ipt1..ipt2) / fpntbu (ipt2..ipt3) /
        // fpntbv (ipt3..iszwr) zones.
        {
            let (w_uhermt, w_rest) = wrkar.split_at_mut(ipt1 as usize);
            let (w_vhermt, w_rest2) = w_rest.split_at_mut(ipt2 as usize - ipt1 as usize);
            let (w_fpntbu, w_fpntbv) =
                w_rest2.split_at_mut(ipt3 as usize - ipt2 as usize);
            mma2cd1_(
                ndimen,
                nbpntu,
                urootl,
                nbpntv,
                vrootl,
                iordru,
                iordrv,
                contr1,
                contr2,
                contr3,
                contr4,
                w_fpntbu,
                w_fpntbv,
                w_uhermt,
                w_vhermt,
                sosotb,
                soditb,
                disotb,
                diditb,
            );
        }
    }

    // ------------------------------ The End -------------------------------
    // OCCT L9999 (iofwr != 0 after a successful allocation: release block).
    {
        let mut isize_ = 0;
        let mut iofst: isize = 0;
        let mut t: Vec<f64> = Vec::new();
        let mut iunit = c__8;
        ier = 0;
        sys_base::mcrdelt_(&mut iunit, &mut isize_, &mut t, &mut iofst, &mut ier);
    }
    if ier > 0 {
        *iercod = 13;
    }
    sys_base::maermsg_("MMA2CDI", iercod);
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA2CDI");
    }
}
