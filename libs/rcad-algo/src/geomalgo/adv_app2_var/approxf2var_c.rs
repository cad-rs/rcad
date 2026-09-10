//! OCCT AdvApp2Var_ApproxF2var (AdvApp2Var_ApproxF2var.cxx) - part C: the
//! iso-discretization engines (mmjacpt_, mma1cdi_, mma1cnt_, mma1fdi_,
//! mma1fer_, mma1jak_, mma1noc_, mma1nop_) plus the
//! [`EvaluatorFunc2Var`] port of AdvApp2Var_EvaluatorFunc2Var.hxx.
//!
//! Encoding note: as in parts A/B, the OCCT Fortran 1-based index formulas
//! are kept on 0-based slices by subtracting the OCCT virtual origin inside
//! the flat-index expression; `&arr[k]` (1-based k) passes as
//! `&arr[(k - 1) as usize..]` (or the raw sub-slice when the base was
//! already offset-adjusted).

use super::approxf2var_b::mmapptt_;
use super::approxf2var_b::mmmapcoe_;
use super::approxf2var_b::mmaperm_;
use super::approxf2var_a::mma1her_;
use super::math_base::pow__di;
use super::math_base_b::mmmpocur_;
use super::math_base_b::mmtrpjj_;
use super::math_base_b::mmaperx_;
use super::sys_base;

// ---------------------------------------------------------------------------
// AdvApp2Var_EvaluatorFunc2Var (AdvApp2Var_EvaluatorFunc2Var.hxx L29-73)
// ---------------------------------------------------------------------------

/// OCCT AdvApp2Var_EvaluatorFunc2Var - the virtual evaluator served to the
/// discretization engines.  The raw-buffer calling convention of
/// `Evaluate` is preserved (AdvApp2Var_EvaluatorFunc2Var.hxx L45-57; the
/// Result layout is Result[Dimension, N]).
pub trait EvaluatorFunc2Var {
    /// OCCT Evaluate (the Evaluate method; inputs are read-only, Result and
    /// ErrorCode are written).
    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        dimension: &i32,
        u_start_end: &[f64],
        v_start_end: &[f64],
        favor_iso: &i32,
        const_param: &f64,
        nb_params: &i32,
        parameters: &[f64],
        u_order: &i32,
        v_order: &i32,
        result: &mut [f64],
        error_code: &mut i32,
    );
}

// ---------------------------------------------------------------------------
// mmjacpt_ (AdvApp2Var_ApproxF2var.cxx L8645-8792)
// ---------------------------------------------------------------------------

/// OCCT mmjacpt_ (AdvApp2Var_ApproxF2var.cxx L8645-8792) - passage from
/// canonical to Jacobi base for a "square" in a space of arbitrary
/// dimension.  All buffers are raw (the OCCT bodies address them through
/// offset-adjusted pointers; the rcad slices are the physical buffers):
/// `ptclgd` / `ptccan` hold (ncoefu, ncoefv, ndimen) columns, `ptcaux`
/// holds 2 * ncoefu * ncoefv * ndimen doubles of workspace.
pub fn mmjacpt_(
    ndimen: &i32,
    ncoefu: &i32,
    ncoefv: &i32,
    iordru: &i32,
    iordrv: &i32,
    ptclgd: &[f64],
    ptcaux: &mut [f64],
    ptccan: &mut [f64],
) {
    // OCCT: ptccan_dim1 = *ncoefu; ptccan_dim2 = *ncoefv;
    //       ptccan_offset = ptccan_dim1 * (ptccan_dim2 + 1) + 1.
    let ptccan_dim1 = *ncoefu;
    let ptccan_dim2 = *ncoefv;
    // OCCT: ptcaux_dim1 = *ncoefv; ptcaux_dim2 = *ncoefu; ptcaux_dim3 = *ndimen.
    let ptcaux_dim1 = *ncoefv;
    let ptcaux_dim2 = *ncoefu;
    let ptcaux_dim3 = *ndimen;

    let pc = |ii: i32, jj: i32, nd: i32| -> usize {
        // OCCT: ptccan[ii + (jj + nd * ptccan_dim2) * ptccan_dim1] - offset
        //   = (ii - 1) + ((jj - 1) + (nd - 1) * dim2) * dim1.
        ((ii - 1) + ((jj - 1) + (nd - 1) * ptccan_dim2) * ptccan_dim1) as usize
    };
    let pa = |jj: i32, ii: i32, ndblk: i32| -> usize {
        // OCCT: ptcaux[jj + (ii + ndblk * ptcaux_dim2) * ptcaux_dim1] - offset
        //   = (jj - 1) + ((ii - 1) + (ndblk - 1) * dim2) * dim1.
        ((jj - 1) + ((ii - 1) + (ndblk - 1) * ptcaux_dim2) * ptcaux_dim1) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMJACPT");
    }

    //   Passage into canonical by u.
    let mut kdim = *ndimen * *ncoefv;
    // OCCT: mmjaccv_(ncoefu, &kdim, iordru, &ptclgd[ptclgd_offset],
    //                &ptcaux[ptcaux_offset], &ptccan[ptccan_offset])
    //   (all three sub-pointers resolve to the raw addresses).
    super::math_base_b::mmjaccv_(ncoefu, &mut kdim, iordru, ptclgd, ptcaux, ptccan);

    //   Swapping of u and v.
    for nd in 1..=*ndimen {
        for jj in 1..=*ncoefv {
            for ii in 1..=*ncoefu {
                // OCCT: ptcaux[jj + (ii + (nd + dim3) * dim2) * dim1]
                //   = ptccan[ii + (jj + nd * dim2) * dim1].
                ptcaux[pa(jj, ii, nd + ptcaux_dim3)] = ptccan[pc(ii, jj, nd)];
            }
        }
    }

    //   Passage into canonical by v.
    kdim = *ndimen * *ncoefu;
    {
        // OCCT: mmjaccv_(ncoefv, &kdim, iordrv,
        //   &ptcaux[((dim3 + 1) * dim2 + 1) * dim1 + 1]        (raw 0),
        //   &ptccan[ptccan_offset]                             (raw 0, scratch),
        //   &ptcaux[(((dim3 << 1) + 1) * dim2 + 1) * dim1 + 1] (raw D)).
        let d = (ptcaux_dim1 * ptcaux_dim2 * ptcaux_dim3) as usize;
        let (ptcaux_in, ptcaux_out) = ptcaux.split_at_mut(d);
        super::math_base_b::mmjaccv_(ncoefv, &mut kdim, iordrv, ptcaux_in, ptccan, ptcaux_out);
    }

    //  Swapping of u and v.
    for nd in 1..=*ndimen {
        for jj in 1..=*ncoefv {
            for ii in 1..=*ncoefu {
                // OCCT: ptccan[ii + (jj + nd * dim2) * dim1]
                //   = ptcaux[jj + (ii + (nd + (dim3 << 1)) * dim2) * dim1].
                ptccan[pc(ii, jj, nd)] = ptcaux[pa(jj, ii, nd + (ptcaux_dim3 << 1))];
            }
        }
    }

    // ---------------------------- THAT'S ALL FOLKS ------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMJACPT");
    }
}

// ---------------------------------------------------------------------------
// mma1cdi_ (AdvApp2Var_ApproxF2var.cxx L233-468)
// ---------------------------------------------------------------------------

/// OCCT mma1cdi_ (AdvApp2Var_ApproxF2var.cxx L233-468) - discretization of
/// the order-IORDRE constraint interpolation polynoms.  `rootlg` is 1-based
/// emulated; the other buffers are raw physical slices.
pub fn mma1cdi_(
    ndimen: &i32,
    nbroot: &i32,
    rootlg: &[f64],
    iordre: &i32,
    contr1: &mut [f64],
    contr2: &mut [f64],
    somtab: &mut [f64],
    diftab: &mut [f64],
    fpntab: &mut [f64],
    hermit: &mut [f64],
    iercod: &mut i32,
) {
    // OCCT: diftab_dim1 = *nbroot / 2 + 1; diftab_offset = diftab_dim1.
    let diftab_dim1 = (nbroot / 2 + 1) as usize;
    let somtab_dim1 = diftab_dim1;
    // OCCT: hermit_dim1 = (*iordre << 1) + 2; hermit_offset = hermit_dim1.
    let hermit_dim1 = ((*iordre << 1) + 2) as usize;
    // OCCT: fpntab_dim1 = *nbroot; fpntab_offset = fpntab_dim1 + 1.
    let fpntab_dim1 = *nbroot;
    let fpntab_off = fpntab_dim1 + 1;
    // OCCT: contr2_dim1 = *ndimen; contr2_offset = contr2_dim1 + 1.
    let contr1_dim1 = *ndimen;
    let contr1_off = contr1_dim1 + 1;
    let contr2_dim1 = *ndimen;
    let contr2_off = contr2_dim1 + 1;

    let st = |ii: i32, nd: i32| -> usize {
        // OCCT: somtab[ii + nd * somtab_dim1] - offset.
        ((ii + nd * somtab_dim1 as i32) - somtab_dim1 as i32) as usize
    };
    let dt = |ii: i32, nd: i32| -> usize {
        ((ii + nd * diftab_dim1 as i32) - diftab_dim1 as i32) as usize
    };
    let ct = |tab_off: i32, dim1: i32, nd: i32, ii: i32| -> usize {
        // OCCT: contrX[nd + ii * contrX_dim1] - offset.
        ((nd + ii * dim1) - tab_off) as usize
    };
    let fp = |kk: i32, ii: i32| -> usize {
        // OCCT: fpntab[kk + ii * fpntab_dim1] - fpntab_offset.
        ((kk + ii * fpntab_dim1) - fpntab_off) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA1CDI");
    }
    *iercod = 0;

    // --- Recuperate 2*(IORDRE+1) coeff of 2*(IORDRE+1) of Hermite polynom ---
    // OCCT: mma1her_(iordre, &hermit[hermit_offset], iercod) (raw address).
    mma1her_(iordre, hermit, iercod);
    if *iercod > 0 {
        // OCCT goto L9100.
        *iercod = 1;
        if ibb >= 3 {
            sys_base::mgsomsg_("MMA1CDI");
        }
        return;
    }

    // ------------------- Discretization of Hermite polynoms ------------
    let mut ncfhe = (*iordre + 1) << 1;
    let mut c__1: i32 = 1;
    for ii in 1..=ncfhe {
        for kk in 1..=*nbroot {
            // OCCT: mmmpocur_(&ncfhe, &c__1, &ncfhe, &hermit[ii * hermit_dim1],
            //                 &rootlg[kk], &fpntab[kk + ii * fpntab_dim1]);
            // hermit sub-pointer raw = (ii - 1) * hermit_dim1;
            // rootlg 1-based emulated; fpntab sub-pointer raw index
            //   kk + ii * dim1 - (dim1 + 1).  OCCT aliases the same int*
            //   twice; the callee only reads NCOFMX/NDEG, so the rcad
            //   encoding splits the alias into two locals.
            let mut ncfhe_a = ncfhe;
            let mut ncfhe_b = ncfhe;
            let mut tparam = rootlg[(kk - 1) as usize];
            mmmpocur_(
                &mut ncfhe_a,
                &mut c__1,
                &mut ncfhe_b,
                &hermit[((ii - 1) * hermit_dim1 as i32) as usize..],
                &mut tparam,
                &mut fpntab[fp(kk, ii)..],
            );
        }
    }

    // ---- Discretizations of boundary polynoms are taken ----
    let nroo2 = *nbroot / 2;
    for nd in 1..=*ndimen {
        for ii in 1..=*iordre + 1 {
            let bid1 = contr1[ct(contr1_off, contr1_dim1, nd, ii)];
            let bid2 = contr2[ct(contr2_off, contr2_dim1, nd, ii)];
            let mut bid3;
            for kk in 1..=nroo2 {
                let kkm = nroo2 - kk + 1;
                bid3 = bid1 * fpntab[fp(kkm, (ii << 1) - 1)] + bid2 * fpntab[fp(kkm, ii << 1)];
                somtab[st(kk, nd)] -= bid3;
                diftab[dt(kk, nd)] += bid3;
            }
            for kk in 1..=nroo2 {
                let kkp = (*nbroot + 1) / 2 + kk;
                bid3 = bid1 * fpntab[fp(kkp, (ii << 1) - 1)] + bid2 * fpntab[fp(kkp, ii << 1)];
                somtab[st(kk, nd)] -= bid3;
                diftab[dt(kk, nd)] -= bid3;
            }
        }
    }

    // ------------ Cas when discretization is done on the roots of a -------
    // ---------- Legendre polynom of uneven degree, 0 is root --------
    if *nbroot % 2 == 1 {
        for nd in 1..=*ndimen {
            let mut bid3 = 0.;
            for ii in 1..=*iordre + 1 {
                bid3 = fpntab[fp(nroo2 + 1, (ii << 1) - 1)] * contr1[ct(contr1_off, contr1_dim1, nd, ii)]
                    + fpntab[fp(nroo2 + 1, ii << 1)] * contr2[ct(contr2_off, contr2_dim1, nd, ii)];
            }
            // OCCT: somtab[nd * somtab_dim1] - offset (ii = 0 slot).
            somtab[st(0, nd)] -= bid3;
            diftab[dt(0, nd)] -= bid3;
        }
    }

    // ------------------------------ The End -------------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA1CDI");
    }
}

// ---------------------------------------------------------------------------
// mma1cnt_ (AdvApp2Var_ApproxF2var.cxx L469-595)
// ---------------------------------------------------------------------------

/// OCCT mma1cnt_ (AdvApp2Var_ApproxF2var.cxx L469-595) - add constraint to
/// polynom.
pub fn mma1cnt_(
    ndimen: &i32,
    iordre: &i32,
    contr1: &[f64],
    contr2: &[f64],
    hermit: &[f64],
    ndgjac: &i32,
    crvjac: &mut [f64],
) {
    // OCCT: hermit_dim1 = (*iordre << 1) + 2; hermit_offset = hermit_dim1.
    let hermit_dim1 = ((*iordre << 1) + 2) as usize;
    // OCCT: contrX_dim1 = *ndimen; contrX_offset = contrX_dim1 + 1.
    let contr1_dim1 = *ndimen;
    let contr1_off = contr1_dim1 + 1;
    let contr2_dim1 = *ndimen;
    let contr2_off = contr2_dim1 + 1;
    // OCCT: crvjac_dim1 = *ndgjac + 1; crvjac_offset = crvjac_dim1.
    let crvjac_dim1 = (*ndgjac + 1) as usize;

    let ct = |tab_off: i32, dim1: i32, nd: i32, jj: i32| -> usize {
        ((nd + jj * dim1) - tab_off) as usize
    };
    let hm = |ii: i32, col: i32| -> usize {
        // OCCT: hermit[ii + col * hermit_dim1] - offset.
        ((ii + col * hermit_dim1 as i32) - hermit_dim1 as i32) as usize
    };
    let cj = |ii: i32, nd: i32| -> usize {
        // OCCT: crvjac[ii + nd * crvjac_dim1] - offset.
        ((ii + nd * crvjac_dim1 as i32) - crvjac_dim1 as i32) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA1CNT");
    }

    for nd in 1..=*ndimen {
        for ii in 0..=(*iordre << 1) + 1 {
            let mut bid = 0.;
            for jj in 1..=*iordre + 1 {
                bid = bid
                    + contr1[ct(contr1_off, contr1_dim1, nd, jj)] * hermit[hm(ii, (jj << 1) - 1)]
                    + contr2[ct(contr2_off, contr2_dim1, nd, jj)] * hermit[hm(ii, jj << 1)];
            }
            crvjac[cj(ii, nd)] = bid;
        }
    }

    if ibb >= 3 {
        sys_base::mgsomsg_("MMA1CNT");
    }
}

// ---------------------------------------------------------------------------
// mma1fdi_ (AdvApp2Var_ApproxF2var.cxx L596-984)
// ---------------------------------------------------------------------------

/// OCCT mma1fdi_ (AdvApp2Var_ApproxF2var.cxx L596-984) - discretization of a
/// non-polynomial function F(U,V) or of its derivative with fixed
/// isoparameter.  `uvfonc` is a raw 4-entry buffer [u0, u1, v0, v1];
/// `ttable` is a raw (nbroot + 2)-entry buffer; `fpntab` is a raw
/// (ndimen, nbp) column buffer of ndimen * nbp doubles; somtab/diftab/
/// contr1/contr2 are raw physical slices.
pub fn mma1fdi_(
    ndimen: &i32,
    uvfonc: &[f64],
    foncnp: &dyn EvaluatorFunc2Var,
    isofav: &i32,
    tconst: &mut f64,
    nbroot: &i32,
    ttable: &mut [f64],
    iordre: &i32,
    ideriv: &i32,
    fpntab: &mut [f64],
    somtab: &mut [f64],
    diftab: &mut [f64],
    contr1: &mut [f64],
    contr2: &mut [f64],
    iercod: &mut i32,
) {
    // OCCT: diftab_dim1 = *nbroot / 2 + 1; diftab_offset = diftab_dim1.
    let diftab_dim1 = (nbroot / 2 + 1) as usize;
    let somtab_dim1 = diftab_dim1;
    // OCCT: fpntab_dim1 = *ndimen (--fpntab).
    let fpntab_dim1 = *ndimen;
    // OCCT: contrX_dim1 = *ndimen; contrX_offset = contrX_dim1 + 1.
    let contr1_dim1 = *ndimen;
    let contr1_off = contr1_dim1 + 1;
    let contr2_dim1 = *ndimen;
    let contr2_off = contr2_dim1 + 1;

    let st = |ii: i32, nd: i32| -> usize {
        ((ii + nd * somtab_dim1 as i32) - somtab_dim1 as i32) as usize
    };
    let dt = |ii: i32, nd: i32| -> usize {
        ((ii + nd * diftab_dim1 as i32) - diftab_dim1 as i32) as usize
    };
    let fp = |nd: i32, col: i32| -> usize {
        // OCCT: fpntab[nd + col * fpntab_dim1] - 1.
        ((nd + col * fpntab_dim1) - 1) as usize
    };
    let ct = |tab_off: i32, dim1: i32, nd: i32, ii: i32| -> usize {
        ((nd + ii * dim1) - tab_off) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA1FDI");
    }
    *iercod = 0;

    // --------------- Definition of the nb of points to calculate ----------
    // --> If constraints, the limits are also taken
    let (ideb, ifin) = if *iordre >= 0 {
        (0, *nbroot + 1)
    // --> Otherwise, only Legendre roots (reframed) are used.
    } else {
        (1, *nbroot)
    };
    // --> Nb of point to calculate.
    let mut nbp = ifin - ideb + 1;
    let nroo2 = *nbroot / 2;

    // --------------- Determination of the order of global derivation ------
    // --> ISOFAV takes only values 1 or 2.
    //    if Iso-U, derive by U of order IDERIV
    let mut ideru;
    let mut iderv;
    let renor;
    if *isofav == 1 {
        ideru = *ideriv;
        iderv = 0;
        let d__1 = (uvfonc[1] - uvfonc[0]) / 2.;
        renor = pow__di(&d__1, ideriv);
    //    if Iso-V, derive by V of order IDERIV
    } else {
        ideru = 0;
        iderv = *ideriv;
        let d__1 = (uvfonc[3] - uvfonc[2]) / 2.;
        renor = pow__di(&d__1, ideriv);
    }

    // ----------- Discretization on roots of the  ---------------
    // ---------------------- Legendre polynom of degree NBROOT -------------
    // OCCT: foncnp.Evaluate(ndimen, &uvfonc[0], &uvfonc[2], isofav, tconst,
    //   &nbp, &ttable[ideb], &ideru, &iderv,
    //   &fpntab[ideb * fpntab_dim1 + 1], iercod).
    foncnp.evaluate(
        ndimen,
        &uvfonc[0..2],
        &uvfonc[2..4],
        isofav,
        tconst,
        &nbp,
        &ttable[ideb as usize..],
        &ideru,
        &iderv,
        &mut fpntab[(ideb * fpntab_dim1) as usize..],
        iercod,
    );
    if *iercod > 0 {
        // OCCT goto L9999.
        finish_mma1fdi_(ibb, iercod);
        return;
    }
    for nd in 1..=*ndimen {
        for ii in 1..=nroo2 {
            let iip = (*nbroot + 1) / 2 + ii;
            let iim = nroo2 - ii + 1;
            let bid1 = fpntab[fp(nd, iim)];
            let bid2 = fpntab[fp(nd, iip)];
            somtab[st(ii, nd)] = renor * (bid2 + bid1);
            diftab[dt(ii, nd)] = renor * (bid2 - bid1);
        }
    }

    // ------------ Case when discretisation is done on roots of a ----
    // ---------- Legendre polynom of uneven degree, 0 is root --------
    if *nbroot % 2 == 1 {
        for nd in 1..=*ndimen {
            somtab[st(0, nd)] = renor * fpntab[fp(nd, nroo2 + 1)];
            diftab[dt(0, nd)] = renor * fpntab[fp(nd, nroo2 + 1)];
        }
    } else {
        for nd in 1..=*ndimen {
            somtab[st(0, nd)] = 0.;
            diftab[dt(0, nd)] = 0.;
        }
    }

    // --------------------- Take into account constraints ----------------
    if *iordre >= 0 {
        // --> Recover already calculated extremities.
        for nd in 1..=*ndimen {
            contr1[ct(contr1_off, contr1_dim1, nd, 1)] = renor * fpntab[fp(nd, 0)];
            contr2[ct(contr2_off, contr2_dim1, nd, 1)] =
                renor * fpntab[fp(nd, *nbroot + 1)];
        }
        // --> Nb of points to calculate/call to FONCNP
        nbp = 1;
        let bid1;
        //    If Iso-U, derive by V till order IORDRE
        if *isofav == 1 {
            // --> Factor of normalisation 1st derivative.
            bid1 = (uvfonc[3] - uvfonc[2]) / 2.;
            for iderv_i in 1..=*iordre {
                iderv = iderv_i;
                // OCCT: Evaluate(..., ttable, &ideru, &iderv,
                //   &contr1[(iderv + 1) * contr1_dim1 + 1], iercod).
                foncnp.evaluate(
                    ndimen,
                    &uvfonc[0..2],
                    &uvfonc[2..4],
                    isofav,
                    tconst,
                    &nbp,
                    ttable,
                    &ideru,
                    &iderv,
                    &mut contr1[ct(contr1_off, contr1_dim1, 1, iderv + 1)..],
                    iercod,
                );
                if *iercod > 0 {
                    finish_mma1fdi_(ibb, iercod);
                    return;
                }
            }
            for iderv_i in 1..=*iordre {
                iderv = iderv_i;
                // OCCT: Evaluate(..., &ttable[*nbroot + 1], &ideru, &iderv,
                //   &contr2[(iderv + 1) * contr2_dim1 + 1], iercod).
                foncnp.evaluate(
                    ndimen,
                    &uvfonc[0..2],
                    &uvfonc[2..4],
                    isofav,
                    tconst,
                    &nbp,
                    &ttable[(*nbroot + 1) as usize..],
                    &ideru,
                    &iderv,
                    &mut contr2[ct(contr2_off, contr2_dim1, 1, iderv + 1)..],
                    iercod,
                );
                if *iercod > 0 {
                    finish_mma1fdi_(ibb, iercod);
                    return;
                }
            }
        //    If Iso-V, derive by U till order IORDRE
        } else {
            // --> Factor of normalization  1st derivative.
            bid1 = (uvfonc[1] - uvfonc[0]) / 2.;
            for ideru_i in 1..=*iordre {
                ideru = ideru_i;
                // OCCT: Evaluate(..., ttable, &ideru, &iderv,
                //   &contr1[(ideru + 1) * contr1_dim1 + 1], iercod).
                foncnp.evaluate(
                    ndimen,
                    &uvfonc[0..2],
                    &uvfonc[2..4],
                    isofav,
                    tconst,
                    &nbp,
                    ttable,
                    &ideru,
                    &iderv,
                    &mut contr1[ct(contr1_off, contr1_dim1, 1, ideru + 1)..],
                    iercod,
                );
                if *iercod > 0 {
                    finish_mma1fdi_(ibb, iercod);
                    return;
                }
            }
            for ideru_i in 1..=*iordre {
                ideru = ideru_i;
                // OCCT: Evaluate(..., &ttable[*nbroot + 1], &ideru, &iderv,
                //   &contr2[(ideru + 1) * contr2_dim1 + 1], iercod).
                foncnp.evaluate(
                    ndimen,
                    &uvfonc[0..2],
                    &uvfonc[2..4],
                    isofav,
                    tconst,
                    &nbp,
                    &ttable[(*nbroot + 1) as usize..],
                    &ideru,
                    &iderv,
                    &mut contr2[ct(contr2_off, contr2_dim1, 1, ideru + 1)..],
                    iercod,
                );
                if *iercod > 0 {
                    finish_mma1fdi_(ibb, iercod);
                    return;
                }
            }        }

        // ------------------------- Normalization of derivatives ----------
        // (The function is redefined on (-1,1)*(-1,1))
        let mut bid2 = renor;
        for ii in 1..=*iordre {
            bid2 = bid1 * bid2;
            for nd in 1..=*ndimen {
                contr1[ct(contr1_off, contr1_dim1, nd, ii + 1)] *= bid2;
                contr2[ct(contr2_off, contr2_dim1, nd, ii + 1)] *= bid2;
            }
        }
    }

    // ------------------------------ The end -------------------------------
    finish_mma1fdi_(ibb, iercod);
}

/// OCCT mma1fdi_ L9999 tail (iercod += 100 + maermsg_ + mgsomsg_).
fn finish_mma1fdi_(ibb: i32, iercod: &mut i32) {
    if *iercod > 0 {
        *iercod += 100;
        sys_base::maermsg_("MMA1FDI", iercod);
    }
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA1FDI");
    }
}

// ---------------------------------------------------------------------------
// mma1fer_ (AdvApp2Var_ApproxF2var.cxx L985-1203)
// ---------------------------------------------------------------------------

/// OCCT mma1fer_ (AdvApp2Var_ApproxF2var.cxx L985-1203) - calculate the
/// degree and the errors of approximation of a border.  `ndimse` /
/// `epsapr` / `errmax` / `errmoy` / `ycvmax` are 1-based emulated;
/// `crvjac` is the raw physical slice.
pub fn mma1fer_(
    ndimen: &i32,
    nbsesp: &i32,
    ndimse: &[i32],
    iordre: &mut i32,
    ndgjac: &i32,
    crvjac: &mut [f64],
    ncflim: &mut i32,
    epsapr: &[f64],
    ycvmax: &mut [f64],
    errmax: &mut [f64],
    errmoy: &mut [f64],
    ncoeff: &mut i32,
    iercod: &mut i32,
) {
    let _ = ndimen;
    // OCCT: crvjac_dim1 = *ndgjac + 1; crvjac_offset = crvjac_dim1.
    let crvjac_dim1 = (*ndgjac + 1) as usize;

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA1FER");
    }
    *iercod = 0;
    let mut idim = 1;
    *ncoeff = 0;
    let mut ncfja = *ndgjac + 1;
    let mut ncfnw = 0;

    // ------------ Calculate the degree of the curve and of the Max error --
    // -------------- of approximation for all sub-spaces --------
    for ii in 1..=*nbsesp {
        let mut ndses = ndimse[(ii - 1) as usize];
        // OCCT aliases the same int* for NCOFMX/NCOEFF (and passes it to
        // several callees); those inputs are read-only in every callee, so
        // the rcad encoding splits the alias into locals.
        let mut ncfja_a = ncfja;
        let mut ncfja_b = ncfja;

        // ------------ cutting of coeff. and calculation of Max error -------
        // OCCT: mmtrpjj_(&ncfja, &ndses, &ncfja, &epsapr[ii], iordre,
        //   &crvjac[idim * crvjac_dim1], &ycvmax[1], &errmax[ii], &ncfnw).
        let mut epsapr_ii = epsapr[(ii - 1) as usize];
        let mut errmax_ii = errmax[(ii - 1) as usize];
        mmtrpjj_(
            &mut ncfja_a,
            &mut ndses,
            &mut ncfja_b,
            &mut epsapr_ii,
            iordre,
            &crvjac[((idim - 1) * crvjac_dim1 as i32) as usize..],
            ycvmax,
            &mut errmax_ii,
            &mut ncfnw,
        );
        errmax[(ii - 1) as usize] = errmax_ii;

        // ------------- If precision OK, calculate the average error -------
        if ncfnw <= *ncflim {
            // OCCT: mmaperm_(&ncfja, &ndses, &ncfja, iordre,
            //   &crvjac[idim * crvjac_dim1], &ncfnw, &errmoy[ii]).
            let mut errmoy_ii = errmoy[(ii - 1) as usize];
            mmaperm_(
                &mut ncfja_a,
                &mut ndses,
                &mut ncfja_b,
                iordre,
                &crvjac[((idim - 1) * crvjac_dim1 as i32) as usize..],
                &mut ncfnw,
                &mut errmoy_ii,
            );
            errmoy[(ii - 1) as usize] = errmoy_ii;
            *ncoeff = ncfnw.max(*ncoeff);

            // ------------- Set the declined coefficients to 0.D0 ----------
            let nbr0 = *ncflim - ncfnw;
            if nbr0 > 0 {
                for kk in 1..=ndses {
                    // OCCT: mvriraz_(&nbr0,
                    //   &crvjac[ncfnw + (idim + kk - 1) * crvjac_dim1]).
                    let raw =
                        ((ncfnw - 1) + ((idim + kk - 1) - 1) * crvjac_dim1 as i32) as usize;
                    sys_base::mvriraz_(nbr0, &mut crvjac[raw..]);
                }
            }
        } else {
            // ------------------- If required precision can't be reached----
            *iercod = -1;

            // ------------------------- calculate the Max error ------------
            // OCCT: mmaperx_(&ncfja, &ndses, &ncfja, iordre,
            //   &crvjac[idim * crvjac_dim1], ncflim, &ycvmax[1],
            //   &errmax[ii], &ier).
            let mut ier = 0;
            let mut errmax_ii = errmax[(ii - 1) as usize];
            mmaperx_(
                &mut ncfja_a,
                &mut ndses,
                &mut ncfja_b,
                iordre,
                &crvjac[((idim - 1) * crvjac_dim1 as i32) as usize..],
                ncflim,
                ycvmax,
                &mut errmax_ii,
                &mut ier,
            );
            errmax[(ii - 1) as usize] = errmax_ii;
            if ier > 0 {
                // OCCT goto L9100.
                *iercod = 1;
                // OCCT goto L9999 (skips idim += ndses).
                break;
            }

            // -------------------- nb of coeff to be returned -------------
            *ncoeff = *ncflim;

            // ------------------- and calculation of the average error ----
            let mut errmoy_ii = errmoy[(ii - 1) as usize];
            mmaperm_(
                &mut ncfja_a,
                &mut ndses,
                &mut ncfja_b,
                iordre,
                &crvjac[((idim - 1) * crvjac_dim1 as i32) as usize..],
                ncflim,
                &mut errmoy_ii,
            );
            errmoy[(ii - 1) as usize] = errmoy_ii;
        }
        idim += ndses;
    }

    // ------------------------------ The end -------------------------------
    // OCCT L9999:
    if *iercod != 0 {
        sys_base::maermsg_("MMA1FER", iercod);
    }
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA1FER");
    }
}

// ---------------------------------------------------------------------------
// mma1jak_ (AdvApp2Var_ApproxF2var.cxx L1368-1489)
// ---------------------------------------------------------------------------

/// OCCT mma1jak_ (AdvApp2Var_ApproxF2var.cxx L1368-1489) - calculate the
/// curve of approximation of a non-polynomial function in the base of
/// Jacobi.  somtab / diftab / cgauss / crvjac are raw physical slices.
pub fn mma1jak_(
    ndimen: &mut i32,
    nbroot: &mut i32,
    iordre: &mut i32,
    ndgjac: &mut i32,
    somtab: &[f64],
    diftab: &[f64],
    cgauss: &mut [f64],
    crvjac: &mut [f64],
    iercod: &mut i32,
) {
    let ibb = sys_base::mnfndeb_();
    if ibb >= 2 {
        sys_base::mgenmsg_("MMA1JAK");
    }
    *iercod = 0;

    // ----------------- Recover coeffs of integration by Gauss -----------
    // OCCT: mmapptt_(ndgjac, nbroot, iordre, cgauss, iercod).
    mmapptt_(ndgjac, nbroot, iordre, cgauss, iercod);
    if *iercod > 0 {
        *iercod = 33;
        // OCCT goto L9999.
        if *iercod != 0 {
            sys_base::maermsg_("MMA1JAK", iercod);
        }
        if ibb >= 2 {
            sys_base::mgsomsg_("MMA1JAK");
        }
        return;
    }

    // --------------- Calculate the curve in the base of Jacobi -----------
    // OCCT: mmmapcoe_(ndimen, ndgjac, iordre, nbroot,
    //   &somtab[somtab_offset], &diftab[diftab_offset], cgauss,
    //   &crvjac[crvjac_offset]) (all sub-pointers raw addresses).
    mmmapcoe_(ndimen, ndgjac, iordre, nbroot, somtab, diftab, cgauss, crvjac);

    // ------------------------------ The End -------------------------------
    if *iercod != 0 {
        sys_base::maermsg_("MMA1JAK", iercod);
    }
    if ibb >= 2 {
        sys_base::mgsomsg_("MMA1JAK");
    }
}

// ---------------------------------------------------------------------------
// mma1noc_ (AdvApp2Var_ApproxF2var.cxx L1490-1626)
// ---------------------------------------------------------------------------

/// OCCT mma1noc_ (AdvApp2Var_ApproxF2var.cxx L1490-1626) - normalization of
/// derivative constraints defined on DFUVIN onto block DUVOUT.  `dfuvin` /
/// `duvout` are raw 4-entry buffers [u0, u1, v0, v1]; `cntrin` / `cntout`
/// are 1-based emulated.
pub fn mma1noc_(
    dfuvin: &[f64],
    ndimen: &i32,
    iordre: &i32,
    cntrin: &[f64],
    duvout: &[f64],
    isofav: &i32,
    ideriv: &i32,
    cntout: &mut [f64],
) {
    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA1NOC");
    }

    // --------------- Determination of coefficients of normalization -------
    let rider;
    let riord;
    if *isofav == 1 {
        // OCCT: d__1 = (dfuvin[4] - dfuvin[3]) / (duvout[4] - duvout[3]).
        let d__1 = (dfuvin[1] - dfuvin[0]) / (duvout[1] - duvout[0]);
        rider = pow__di(&d__1, ideriv);
        // OCCT: d__1 = (dfuvin[6] - dfuvin[5]) / (duvout[6] - duvout[5]).
        let d__1 = (dfuvin[3] - dfuvin[2]) / (duvout[3] - duvout[2]);
        riord = pow__di(&d__1, iordre);
    } else {
        let d__1 = (dfuvin[3] - dfuvin[2]) / (duvout[3] - duvout[2]);
        rider = pow__di(&d__1, ideriv);
        let d__1 = (dfuvin[1] - dfuvin[0]) / (duvout[1] - duvout[0]);
        riord = pow__di(&d__1, iordre);
    }

    // ------------- Renormalization of the vector of constraint ------------
    let bid = rider * riord;
    for nd in 1..=*ndimen {
        // OCCT: cntout[nd] = bid * cntrin[nd] (1-based).
        cntout[(nd - 1) as usize] = bid * cntrin[(nd - 1) as usize];
    }

    // ------------------------------ The end -------------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA1NOC");
    }
}

// ---------------------------------------------------------------------------
// mma1nop_ (AdvApp2Var_ApproxF2var.cxx L1627-1736)
// ---------------------------------------------------------------------------

/// OCCT mma1nop_ (AdvApp2Var_ApproxF2var.cxx L1627-1736) - normalization of
/// the parameters of an iso from the parametric block and the parameters on
/// (-1,1).  `rootlg` is 1-based emulated; `uvfonc` / `ttable` are raw
/// buffers (ttable has nbroot + 2 entries).
pub fn mma1nop_(
    nbroot: &i32,
    rootlg: &[f64],
    uvfonc: &[f64],
    isofav: &i32,
    ttable: &mut [f64],
    iercod: &mut i32,
) {
    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA1NOP");
    }

    // OCCT: alinu = (uvfonc[4] - uvfonc[3]) / 2. (adjusted -= 3 = raw [0..1]).
    let alinu = (uvfonc[1] - uvfonc[0]) / 2.;
    let blinu = (uvfonc[1] + uvfonc[0]) / 2.;
    // OCCT: alinv = (uvfonc[6] - uvfonc[5]) / 2. (raw [2..3]).
    let alinv = (uvfonc[3] - uvfonc[2]) / 2.;
    let blinv = (uvfonc[3] + uvfonc[2]) / 2.;

    if *isofav == 1 {
        // OCCT: ttable[0] = uvfonc[5] (raw [2]).
        ttable[0] = uvfonc[2];
        for ii in 1..=*nbroot {
            // OCCT: ttable[ii] = alinv * rootlg[ii] + blinv.
            ttable[ii as usize] = alinv * rootlg[(ii - 1) as usize] + blinv;
        }
        ttable[(*nbroot + 1) as usize] = uvfonc[3];
    } else if *isofav == 2 {
        ttable[0] = uvfonc[0];
        for ii in 1..=*nbroot {
            ttable[ii as usize] = alinu * rootlg[(ii - 1) as usize] + blinu;
        }
        ttable[(*nbroot + 1) as usize] = uvfonc[1];
    } else {
        // OCCT goto L9100.
        *iercod = 1;
        if *iercod != 0 {
            sys_base::maermsg_("MMA1NOP", iercod);
        }
        if ibb >= 3 {
            sys_base::mgsomsg_("MMA1NOP");
        }
        return;
    }

    // ------------------------------ THE END -------------------------------
    if *iercod != 0 {
        sys_base::maermsg_("MMA1NOP", iercod);
    }
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA1NOP");
    }
}
