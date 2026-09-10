//! OCCT AdvApp2Var_ApproxF2var (AdvApp2Var_ApproxF2var.cxx) - part A: the
//! Hermit/edge-constraint engines consumed by AdvApp2Var_Patch::AddConstraints.
//!
//! Encoding note: the OCCT engines use Fortran 1-based storage with
//! `pointer -= offset` adjustments; rcad keeps the OCCT index formulas on
//! 0-based slices by subtracting the OCCT virtual origin inside the
//! flat-index expression (annotated per function as `*_offset`).

use super::sys_base;

/// OCCT AdvApp2Var_ApproxF2var::mma1her_ (AdvApp2Var_ApproxF2var.cxx
/// L1204-1364) - the 2*(IORDRE+1) Hermit polynomials of degree 2*IORDRE+1 on
/// (-1,1).  `hermit` is dimensioned (2*(IORDRE+1), 2*(IORDRE+1)) Fortran-wise;
/// the rcad slice is the physical (1-based emulated) buffer.
pub fn mma1her_(iordre: &i32, hermit: &mut [f64], iercod: &mut i32) {
    // OCCT: hermit_dim1 = (iordre + 1) << 1; hermit_offset = dim1 + 1.
    let hermit_dim1 = ((*iordre + 1) << 1) as usize;
    let hermit_off = hermit_dim1 + 1;
    let h = |virtual_index: usize| -> usize {
        virtual_index - hermit_off
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA1HER");
    }
    *iercod = 0;

    // --- Recover (IORDRE+2) coeff of 2*(IORDRE+1) Hermit polynoms ---
    if *iordre == 0 {
        hermit[h(hermit_dim1 + 1)] = 0.5;
        hermit[h(hermit_dim1 + 2)] = -0.5;

        hermit[h((hermit_dim1 << 1) + 1)] = 0.5;
        hermit[h((hermit_dim1 << 1) + 2)] = 0.5;
    } else if *iordre == 1 {
        hermit[h(hermit_dim1 + 1)] = 0.5;
        hermit[h(hermit_dim1 + 2)] = -0.75;
        hermit[h(hermit_dim1 + 3)] = 0.0;
        hermit[h(hermit_dim1 + 4)] = 0.25;

        hermit[h((hermit_dim1 << 1) + 1)] = 0.5;
        hermit[h((hermit_dim1 << 1) + 2)] = 0.75;
        hermit[h((hermit_dim1 << 1) + 3)] = 0.0;
        hermit[h((hermit_dim1 << 1) + 4)] = -0.25;

        hermit[h(hermit_dim1 * 3 + 1)] = 0.25;
        hermit[h(hermit_dim1 * 3 + 2)] = -0.25;
        hermit[h(hermit_dim1 * 3 + 3)] = -0.25;
        hermit[h(hermit_dim1 * 3 + 4)] = 0.25;

        hermit[h((hermit_dim1 << 2) + 1)] = -0.25;
        hermit[h((hermit_dim1 << 2) + 2)] = -0.25;
        hermit[h((hermit_dim1 << 2) + 3)] = 0.25;
        hermit[h((hermit_dim1 << 2) + 4)] = 0.25;
    } else if *iordre == 2 {
        hermit[h(hermit_dim1 + 1)] = 0.5;
        hermit[h(hermit_dim1 + 2)] = -0.9375;
        hermit[h(hermit_dim1 + 3)] = 0.0;
        hermit[h(hermit_dim1 + 4)] = 0.625;
        hermit[h(hermit_dim1 + 5)] = 0.0;
        hermit[h(hermit_dim1 + 6)] = -0.1875;

        hermit[h((hermit_dim1 << 1) + 1)] = 0.5;
        hermit[h((hermit_dim1 << 1) + 2)] = 0.9375;
        hermit[h((hermit_dim1 << 1) + 3)] = 0.0;
        hermit[h((hermit_dim1 << 1) + 4)] = -0.625;
        hermit[h((hermit_dim1 << 1) + 5)] = 0.0;
        hermit[h((hermit_dim1 << 1) + 6)] = 0.1875;

        hermit[h(hermit_dim1 * 3 + 1)] = 0.3125;
        hermit[h(hermit_dim1 * 3 + 2)] = -0.4375;
        hermit[h(hermit_dim1 * 3 + 3)] = -0.375;
        hermit[h(hermit_dim1 * 3 + 4)] = 0.625;
        hermit[h(hermit_dim1 * 3 + 5)] = 0.0625;
        hermit[h(hermit_dim1 * 3 + 6)] = -0.1875;

        hermit[h((hermit_dim1 << 2) + 1)] = -0.3125;
        hermit[h((hermit_dim1 << 2) + 2)] = -0.4375;
        hermit[h((hermit_dim1 << 2) + 3)] = 0.375;
        hermit[h((hermit_dim1 << 2) + 4)] = 0.625;
        hermit[h((hermit_dim1 << 2) + 5)] = -0.0625;
        hermit[h((hermit_dim1 << 2) + 6)] = -0.1875;

        hermit[h(hermit_dim1 * 5 + 1)] = 0.0625;
        hermit[h(hermit_dim1 * 5 + 2)] = -0.0625;
        hermit[h(hermit_dim1 * 5 + 3)] = -0.125;
        hermit[h(hermit_dim1 * 5 + 4)] = 0.125;
        hermit[h(hermit_dim1 * 5 + 5)] = 0.0625;
        hermit[h(hermit_dim1 * 5 + 6)] = -0.0625;

        hermit[h(hermit_dim1 * 6 + 1)] = 0.0625;
        hermit[h(hermit_dim1 * 6 + 2)] = 0.0625;
        hermit[h(hermit_dim1 * 6 + 3)] = -0.125;
        hermit[h(hermit_dim1 * 6 + 4)] = -0.125;
        hermit[h(hermit_dim1 * 6 + 5)] = 0.0625;
        hermit[h(hermit_dim1 * 6 + 6)] = 0.0625;
    } else {
        *iercod = 1;
    }

    // ------------------------------ The End -------------------------------
    sys_base::maermsg_("MMA1HER", iercod);
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA1HER");
    }
}

/// OCCT AdvApp2Var_ApproxF2var::mma2ac1_ (AdvApp2Var_ApproxF2var.cxx
/// L1738-1904) - subtraction of the angular constraints from PATJAC.
pub fn mma2ac1_(
    ndimen: &i32,
    mxujac: &i32,
    mxvjac: &i32,
    iordru: &i32,
    iordrv: &i32,
    contr1: &[f64],
    contr2: &[f64],
    contr3: &[f64],
    contr4: &[f64],
    uhermt: &[f64],
    vhermt: &[f64],
    patjac: &mut [f64],
) {
    let patjac_dim1 = (*mxujac + 1) as usize;
    let patjac_dim2 = (*mxvjac + 1) as usize;
    // OCCT: patjac_offset = patjac_dim1 * patjac_dim2.
    let patjac_off = patjac_dim1 * patjac_dim2;
    let uhermt_dim1 = ((*iordru << 1) + 2) as usize;
    // OCCT: uhermt_offset = uhermt_dim1.
    let uhermt_off = uhermt_dim1;
    let vhermt_dim1 = ((*iordrv << 1) + 2) as usize;
    // OCCT: vhermt_offset = vhermt_dim1.
    let vhermt_off = vhermt_dim1;
    let contr_dim1 = *ndimen as usize;
    let contr_dim2 = (*iordru + 2) as usize;
    // OCCT: contrX_offset = contrX_dim1 * (contrX_dim2 + 1) + 1.
    let contr_off = contr_dim1 * (contr_dim2 + 1) + 1;

    let pj = |ku: i32, kv: i32, nd: i32| -> usize {
        // OCCT: patjac[ku + (kv + nd * patjac_dim2) * patjac_dim1] - offset.
        ((ku + (kv + nd * patjac_dim2 as i32) * patjac_dim1 as i32) as usize) - patjac_off
    };
    let ct = |c: &[f64], nd: i32, ii: i32, jj: i32| -> usize {
        // OCCT: contrX[nd + (ii + jj * contrX_dim2) * contrX_dim1] - offset.
        ((nd + (ii + jj * contr_dim2 as i32) * contr_dim1 as i32) as usize) - contr_off
    };
    let uh = |ku: i32, term: i32| -> usize {
        // OCCT: uhermt[ku + term * uhermt_dim1] - uhermt_offset.
        ((ku + term * uhermt_dim1 as i32) as usize) - uhermt_off
    };
    let vh = |kv: i32, term: i32| -> usize {
        // OCCT: vhermt[kv + term * vhermt_dim1] - vhermt_offset.
        ((kv + term * vhermt_dim1 as i32) as usize) - vhermt_off
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2AC1");
    }

    // ------------ SUBTRACTION OF ANGULAR CONSTRAINTS -------------------
    let ioru1 = *iordru + 1;
    let iorv1 = *iordrv + 1;
    let ndgu = (*iordru << 1) + 1;
    let ndgv = (*iordrv << 1) + 1;

    for jj in 1..=iorv1 {
        for ii in 1..=ioru1 {
            for nd in 1..=*ndimen {
                let cnt1 = contr1[ct(contr1, nd, ii, jj)];
                let cnt2 = contr2[ct(contr2, nd, ii, jj)];
                let cnt3 = contr3[ct(contr3, nd, ii, jj)];
                let cnt4 = contr4[ct(contr4, nd, ii, jj)];
                for kv in 0..=ndgv {
                    let bidv1 = vhermt[vh(kv, (jj << 1) - 1)];
                    let bidv2 = vhermt[vh(kv, jj << 1)];
                    for ku in 0..=ndgu {
                        let bidu1 = uhermt[uh(ku, (ii << 1) - 1)];
                        let bidu2 = uhermt[uh(ku, ii << 1)];
                        patjac[pj(ku, kv, nd)] = patjac[pj(ku, kv, nd)]
                            - bidu1 * bidv1 * cnt1
                            - bidu2 * bidv1 * cnt2
                            - bidu1 * bidv2 * cnt3
                            - bidu2 * bidv2 * cnt4;
                    }
                }
            }
        }
    }

    // ------------------------------ The end -------------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2AC1");
    }
}

/// OCCT AdvApp2Var_ApproxF2var::mma2ac2_ (AdvApp2Var_ApproxF2var.cxx
/// L1908-2050) - adding the iso-V constraint coefficients (curves by u,
/// Hermit by v) into PATJAC.
pub fn mma2ac2_(
    ndimen: &i32,
    mxujac: &i32,
    mxvjac: &i32,
    iordrv: &i32,
    nclimu: &i32,
    ncfiv1: &[i32],
    crbiv1: &[f64],
    ncfiv2: &[i32],
    crbiv2: &[f64],
    vhermt: &[f64],
    patjac: &mut [f64],
) {
    let patjac_dim1 = (*mxujac + 1) as usize;
    let patjac_dim2 = (*mxvjac + 1) as usize;
    // OCCT: patjac_offset = patjac_dim1 * patjac_dim2.
    let patjac_off = patjac_dim1 * patjac_dim2;
    let vhermt_dim1 = ((*iordrv << 1) + 2) as usize;
    // OCCT: vhermt_offset = vhermt_dim1.
    let vhermt_off = vhermt_dim1;
    // OCCT: --ncfiv1; --ncfiv2 (1-based slices).
    let crb_dim1 = *nclimu as usize;
    let crb_dim2 = *ndimen as usize;
    // OCCT: crbivX_offset = crbivX_dim1 * (crbivX_dim2 + 1).
    let crb_off = crb_dim1 * (crb_dim2 + 1);

    let pj = |kk: i32, jj: i32, nd: i32| -> usize {
        // OCCT: patjac[kk + (jj + nd * patjac_dim2) * patjac_dim1] - offset.
        ((kk + (jj + nd * patjac_dim2 as i32) * patjac_dim1 as i32) as usize) - patjac_off
    };
    let crb = |c: &[f64], kk: i32, nd: i32, ii: i32| -> usize {
        // OCCT: crbivX[kk + (nd + ii * crbivX_dim2) * crbivX_dim1] - offset.
        ((kk + (nd + ii * crb_dim2 as i32) * crb_dim1 as i32) as usize) - crb_off
    };
    let vh = |jj: i32, term: i32| -> usize {
        // OCCT: vhermt[jj + term * vhermt_dim1] - vhermt_offset.
        ((jj + term * vhermt_dim1 as i32) as usize) - vhermt_off
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2AC2");
    }

    // ------------ ADDING of coeff by u of curves, by v of Hermit --------
    for ii in 1..=*iordrv + 1 {
        let ndgv1 = ncfiv1[ii as usize - 1] - 1;
        let ndgv2 = ncfiv2[ii as usize - 1] - 1;
        for nd in 1..=*ndimen {
            for jj in 0..=(*iordrv << 1) + 1 {
                let bid1 = vhermt[vh(jj, (ii << 1) - 1)];
                for kk in 0..=ndgv1 {
                    patjac[pj(kk, jj, nd)] += bid1 * crbiv1[crb(crbiv1, kk, nd, ii)];
                }
                let bid2 = vhermt[vh(jj, ii << 1)];
                for kk in 0..=ndgv2 {
                    patjac[pj(kk, jj, nd)] += bid2 * crbiv2[crb(crbiv2, kk, nd, ii)];
                }
            }
        }
    }

    // ------------------------------ The end -------------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2AC2");
    }
}

/// OCCT AdvApp2Var_ApproxF2var::mma2ac3_ (AdvApp2Var_ApproxF2var.cxx
/// L2054-2201) - adding the iso-U constraint coefficients (curves by v,
/// Hermit by u) into PATJAC.
pub fn mma2ac3_(
    ndimen: &i32,
    mxujac: &i32,
    mxvjac: &i32,
    iordru: &i32,
    nclimv: &i32,
    ncfiu1: &[i32],
    crbiu1: &[f64],
    ncfiu2: &[i32],
    crbiu2: &[f64],
    uhermt: &[f64],
    patjac: &mut [f64],
) {
    let patjac_dim1 = (*mxujac + 1) as usize;
    let patjac_dim2 = (*mxvjac + 1) as usize;
    // OCCT: patjac_offset = patjac_dim1 * patjac_dim2.
    let patjac_off = patjac_dim1 * patjac_dim2;
    let uhermt_dim1 = ((*iordru << 1) + 2) as usize;
    // OCCT: uhermt_offset = uhermt_dim1.
    let uhermt_off = uhermt_dim1;
    // OCCT: --ncfiu1; --ncfiu2 (1-based slices).
    let crb_dim1 = *nclimv as usize;
    let crb_dim2 = *ndimen as usize;
    // OCCT: crbiuX_offset = crbiuX_dim1 * (crbiuX_dim2 + 1).
    let crb_off = crb_dim1 * (crb_dim2 + 1);

    let pj = |kk: i32, jj: i32, nd: i32| -> usize {
        // OCCT: patjac[kk + (jj + nd * patjac_dim2) * patjac_dim1] - offset.
        ((kk + (jj + nd * patjac_dim2 as i32) * patjac_dim1 as i32) as usize) - patjac_off
    };
    let crb = |c: &[f64], jj: i32, nd: i32, ii: i32| -> usize {
        // OCCT: crbiuX[jj + (nd + ii * crbiuX_dim2) * crbiuX_dim1] - offset.
        ((jj + (nd + ii * crb_dim2 as i32) * crb_dim1 as i32) as usize) - crb_off
    };
    let uh = |kk: i32, term: i32| -> usize {
        // OCCT: uhermt[kk + term * uhermt_dim1] - uhermt_offset.
        ((kk + term * uhermt_dim1 as i32) as usize) - uhermt_off
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2AC3");
    }

    // ------------ ADDING of coeff by u of curves, by v of Hermit --------
    for ii in 1..=*iordru + 1 {
        let ndgu1 = ncfiu1[ii as usize - 1] - 1;
        let ndgu2 = ncfiu2[ii as usize - 1] - 1;
        for nd in 1..=*ndimen {
            for jj in 0..=ndgu1 {
                let bid1 = crbiu1[crb(crbiu1, jj, nd, ii)];
                for kk in 0..=(*iordru << 1) + 1 {
                    patjac[pj(kk, jj, nd)] += bid1 * uhermt[uh(kk, (ii << 1) - 1)];
                }
            }
            for jj in 0..=ndgu2 {
                let bid2 = crbiu2[crb(crbiu2, jj, nd, ii)];
                for kk in 0..=(*iordru << 1) + 1 {
                    patjac[pj(kk, jj, nd)] += bid2 * uhermt[uh(kk, ii << 1)];
                }
            }
        }
    }

    // ------------------------------ The end -------------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2AC3");
    }
}
