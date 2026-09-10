//! OCCT AdvApp2Var_ApproxF2var (AdvApp2Var_ApproxF2var.cxx) - part B: the
//! Jacobi/Legendre discretization engines consumed by AdvApp2Var_Context and
//! the coefficient-collection leaves used by the later mma2ce*/mma2fx* set.
//!
//! Encoding note: as in part A, the OCCT Fortran index formulas are kept on
//! 0-based slices with the OCCT virtual origin subtracted inside the
//! flat-index expression; `&arr[k]` (1-based k) passes as
//! `&slice[(k - 1) as usize..]`.

use super::data_mmapgs0::GSL0J0;
use super::data_mmapgs0::GSLXJ0;
use super::data_mmapgs1::GSL0J1;
use super::data_mmapgs1::GSLXJ1;
use super::data_mmapgs2::GSL0J2;
use super::data_mmapgs2::GSLXJ2;
use super::data_mmapgss::GSL0JS;
use super::data_mmapgss::GSLXJS;
use super::math_base::mmrtptt_;
use super::sys_base;

// OCCT mma2jmx_ initialized data (AdvApp2Var_ApproxF2var.cxx L7580-7663).
#[rustfmt::skip]
pub static XMAX2: [f64; 57] = [
    0.9682458365518542212948163499456, 0.986013297183269340427888048593603,
    1.07810420343739860362585159028115, 1.17325804490920057010925920756025,
    1.26476561266905634732910520370741, 1.35169950227289626684434056681946,
    1.43424378958284137759129885012494, 1.51281316274895465689402798226634,
    1.5878364329591908800533936587012, 1.65970112228228167018443636171226,
    1.72874345388622461848433443013543, 1.7952515611463877544077632304216,
    1.85947199025328260370244491818047, 1.92161634324190018916351663207101,
    1.98186713586472025397859895825157, 2.04038269834980146276967984252188,
    2.09730119173852573441223706382076, 2.15274387655763462685970799663412,
    2.20681777186342079455059961912859, 2.25961782459354604684402726624239,
    2.31122868752403808176824020121524, 2.36172618435386566570998793688131,
    2.41117852396114589446497298177554, 2.45964731268663657873849811095449,
    2.50718840313973523778244737914028, 2.55385260994795361951813645784034,
    2.59968631659221867834697883938297, 2.64473199258285846332860663371298,
    2.68902863641518586789566216064557, 2.73261215675199397407027673053895,
    2.77551570192374483822124304745691, 2.8177699459714315371037628127545,
    2.85940333797200948896046563785957, 2.90044232019793636101516293333324,
    2.94091151970640874812265419871976, 2.98083391718088702956696303389061,
    3.02023099621926980436221568258656, 3.05912287574998661724731962377847,
    3.09752842783622025614245706196447, 3.13546538278134559341444834866301,
    3.17295042316122606504398054547289, 3.2099992681699613513775259670214,
    3.24662674946606137764916854570219, 3.28284687953866689817670991319787,
    3.31867291347259485044591136879087, 3.35411740487202127264475726990106,
    3.38919225660177218727305224515862, 3.42390876691942143189170489271753,
    3.45827767149820230182596660024454, 3.49230918177808483937957161007792,
    3.5260130200285724149540352829756, 3.55939845146044235497103883695448,
    3.59247431368364585025958062194665, 3.62524904377393592090180712976368,
    3.65773070318071087226169680450936, 3.68992700068237648299565823810245,
    3.72184531357268220291630708234186,
];
#[rustfmt::skip]
pub static XMAX4: [f64; 55] = [
    1.1092649593311780079813740546678, 1.05299572648705464724876659688996,
    1.0949715351434178709281698645813, 1.15078388379719068145021100764647,
    1.2094863084718701596278219811869, 1.26806623151369531323304177532868,
    1.32549784426476978866302826176202, 1.38142537365039019558329304432581,
    1.43575531950773585146867625840552, 1.48850442653629641402403231015299,
    1.53973611681876234549146350844736, 1.58953193485272191557448229046492,
    1.63797820416306624705258190017418, 1.68515974143594899185621942934906,
    1.73115699602477936547107755854868, 1.77604489805513552087086912113251,
    1.81989256661534438347398400420601, 1.86276344480103110090865609776681,
    1.90471563564740808542244678597105, 1.94580231994751044968731427898046,
    1.98607219357764450634552790950067, 2.02556989246317857340333585562678,
    2.06433638992049685189059517340452, 2.10240936014742726236706004607473,
    2.13982350649113222745523925190532, 2.17661085564771614285379929798896,
    2.21280102016879766322589373557048, 2.2484214321456956597803794333791,
    2.28349755104077956674135810027654, 2.31805304852593774867640120860446,
    2.35210997297725685169643559615022, 2.38568889602346315560143377261814,
    2.41880904328694215730192284109322, 2.45148841120796359750021227795539,
    2.48374387161372199992570528025315, 2.5155912654873773953959098501893,
    2.54704548720896557684101746505398, 2.57812056037881628390134077704127,
    2.60882970619319538196517982945269, 2.63918540521920497868347679257107,
    2.66919945330942891495458446613851, 2.69888301230439621709803756505788,
    2.72824665609081486732853370048, 2.75730041251405791603760003778285,
    2.78605380158311346185098508516203, 2.81451587035387403267676338931454,
    2.84269522483114290814009184272637, 2.87060005919012917988363332454033,
    2.89823818258367657739520912946934, 2.92561704377132528239806135133273,
    2.95274375377994262301217318010209, 2.97962510678256471794289060402033,
    3.00626759936182712291041810228171, 3.03267744830655121818899164295959,
    3.05886060707437081434964933864149,
];
#[rustfmt::skip]
pub static XMAX6: [f64; 53] = [
    1.21091229812484768570102219548814, 1.11626917091567929907256116528817,
    1.1327140810290884106278510474203, 1.1679452722668028753522098022171,
    1.20910611986279066645602153641334, 1.25228283758701572089625983127043,
    1.29591971597287895911380446311508, 1.3393138157481884258308028584917,
    1.3821288728999671920677617491385, 1.42420414683357356104823573391816,
    1.46546895108549501306970087318319, 1.50590085198398789708599726315869,
    1.54550385142820987194251585145013, 1.58429644271680300005206185490937,
    1.62230484071440103826322971668038, 1.65955905239130512405565733793667,
    1.69609056468292429853775667485212, 1.73193098017228915881592458573809,
    1.7671112206990325429863426635397, 1.80166107681586964987277458875667,
    1.83560897003644959204940535551721, 1.86898184653271388435058371983316,
    1.90180515174518670797686768515502, 1.93410285411785808749237200054739,
    1.96589749778987993293150856865539, 1.99721027139062501070081653790635,
    2.02806108474738744005306947877164, 2.05846864831762572089033752595401,
    2.08845055210580131460156962214748, 2.11802334209486194329576724042253,
    2.14720259305166593214642386780469, 2.17600297710595096918495785742803,
    2.20443832785205516555772788192013, 2.2325216999457379530416998244706,
    2.2602654243075083168599953074345, 2.28768115912702794202525264301585,
    2.3147799369092684021274946755348, 2.34157220782483457076721300512406,
    2.36806787963276257263034969490066, 2.39427635443992520016789041085844,
    2.42020656255081863955040620243062, 2.44586699364757383088888037359254,
    2.47126572552427660024678584642791, 2.49641045058324178349347438430311,
    2.52130850028451113942299097584818, 2.54596686772399937214920135190177,
    2.5703922285006754089328998222275, 2.59459096001908861492582631591134,
    2.61856915936049852435394597597773, 2.64233265984385295286445444361827,
    2.66588704638685848486056711408168, 2.68923766976735295746679957665724,
    2.71238965987606292679677228666411,
];

/// OCCT AdvApp2Var_ApproxF2var::mma2jmx_ (AdvApp2Var_ApproxF2var.cxx
/// L7576-7763) - maximums of the Jacobi polynomials times the weight on
/// (-1,1) for orders 0, 4, 6 or Legendre.
pub fn mma2jmx_(ndgjac: &mut i32, iordre: &mut i32, xjacmx: &mut [f64]) {
    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2JMX");
    }

    let numax = *ndgjac - ((*iordre + 1) << 1);
    if *iordre == -1 {
        for ii in 0..=numax {
            let bid = (ii as f64 * 2. + 1.) / 2.;
            xjacmx[ii as usize] = bid.sqrt();
        }
    } else if *iordre == 0 {
        for ii in 0..=numax {
            xjacmx[ii as usize] = XMAX2[ii as usize];
        }
    } else if *iordre == 1 {
        for ii in 0..=numax {
            xjacmx[ii as usize] = XMAX4[ii as usize];
        }
    } else if *iordre == 2 {
        for ii in 0..=numax {
            xjacmx[ii as usize] = XMAX6[ii as usize];
        }
    }

    // ------------------------- The end ------------------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2JMX");
    }
}

/// OCCT AdvApp2Var_ApproxF2var::mma2moy_ (AdvApp2Var_ApproxF2var.cxx
/// L7767-7921) - average error upper bound when only the PATJAC coefficients
/// of degree between 2*(IORDRU+1)..MINDGU / 2*(IORDRV+1)..MINDGV are kept.
pub fn mma2moy_(
    ndgumx: &mut i32,
    ndgvmx: &mut i32,
    ndimen: &mut i32,
    mindgu: &mut i32,
    maxdgu: &mut i32,
    mindgv: &mut i32,
    maxdgv: &mut i32,
    iordru: &mut i32,
    iordrv: &mut i32,
    patjac: &[f64],
    errmoy: &mut f64,
) {
    let patjac_dim1 = (*ndgumx + 1) as usize;
    let patjac_dim2 = (*ndgvmx + 1) as usize;
    // OCCT: patjac_offset = patjac_dim1 * patjac_dim2.
    let patjac_off = (patjac_dim1 * patjac_dim2) as i32;
    let pj = |ii: i32, jj: i32, nd: i32| -> usize {
        // OCCT: patjac[ii + (jj + nd * patjac_dim2) * patjac_dim1] - offset.
        ((ii + (jj + nd * patjac_dim2 as i32) * patjac_dim1 as i32) as usize)
            - patjac_off as usize
    };

    let ldbg = sys_base::mnfndeb_() >= 3;
    if ldbg {
        sys_base::mgenmsg_("MMA2MOY");
    }

    let idebu = (*iordru + 1) << 1;
    let idebv = (*iordrv + 1) << 1;
    let minu = idebu.max(*mindgu);
    let minv = idebv.max(*mindgv);
    let mut bid0 = 0.;
    *errmoy = 0.;

    // ------- upper bound when coeff MINDGU..MAXDGU / MINDGV..MAXDGV removed ----
    for nd in 1..=*ndimen {
        for jj in minv..=*maxdgv {
            for ii in idebu..=*maxdgu {
                let bid1 = patjac[pj(ii, jj, nd)];
                bid0 += bid1 * bid1;
            }
        }
    }

    for nd in 1..=*ndimen {
        for jj in idebv..=minv - 1 {
            for ii in minu..=*maxdgu {
                let bid1 = patjac[pj(ii, jj, nd)];
                bid0 += bid1 * bid1;
            }
        }
    }

    // ----------------------- Calculation of the average error -------------
    bid0 /= 4.;
    *errmoy = bid0.sqrt();

    // ------------------------- The end ------------------------------------
    if ldbg {
        sys_base::mgsomsg_("MMA2MOY");
    }
}

/// OCCT AdvApp2Var_ApproxF2var::mma2roo_ (AdvApp2Var_ApproxF2var.cxx
/// L7925-8024) - roots of Legendre for the discretizations.  `urootl` /
/// `vrootl` carry 1-based OCCT indexing: the OCCT sub-pointer
/// `&arr[k]` maps to `&arr[(k - 1) as usize..]` (the raw address), and the
/// negation / zero writes map their 1-based indices with `- 1`.
pub fn mma2roo_(nbpntu: &mut i32, nbpntv: &mut i32, urootl: &mut [f64], vrootl: &mut [f64]) {
    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA2ROO");
    }

    // ---------------- Return the POSITIVE roots on U ------------------
    // OCCT: mmrtptt_(nbpntu, &urootl[(*nbpntu + 1) / 2 + 1]).
    mmrtptt_(
        nbpntu,
        &mut urootl[((*nbpntu + 1) / 2 + 1 - 1) as usize..],
    );
    for ii in 1..=*nbpntu / 2 {
        urootl[(ii - 1) as usize] = -urootl[(*nbpntu - ii + 1 - 1) as usize];
    }
    if *nbpntu % 2 == 1 {
        urootl[(*nbpntu / 2 + 1 - 1) as usize] = 0.;
    }

    // ---------------- Return the POSITIVE roots on V ------------------
    mmrtptt_(
        nbpntv,
        &mut vrootl[((*nbpntv + 1) / 2 + 1 - 1) as usize..],
    );
    for ii in 1..=*nbpntv / 2 {
        vrootl[(ii - 1) as usize] = -vrootl[(*nbpntv - ii + 1 - 1) as usize];
    }
    if *nbpntv % 2 == 1 {
        vrootl[(*nbpntv / 2 + 1 - 1) as usize] = 0.;
    }

    // ------------------------------ The End -------------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA2ROO");
    }
}

/// OCCT mmmapcoe_ (AdvApp2Var_ApproxF2var.cxx L8028-8194) - least-squares
/// coefficients in the Jacobi base from the Legendre-root discretization.
pub fn mmmapcoe_(
    ndim: &mut i32,
    ndgjac: &mut i32,
    iordre: &mut i32,
    nbpnts: &mut i32,
    somtab: &[f64],
    diftab: &[f64],
    gsstab: &[f64],
    crvjac: &mut [f64],
) {
    let crvjac_dim1 = (*ndgjac + 1) as i32;
    // OCCT: crvjac_offset = crvjac_dim1.
    let crvjac_off = crvjac_dim1;
    let gsstab_dim1 = *nbpnts / 2 + 1;
    let diftab_dim1 = *nbpnts / 2 + 1;
    // OCCT: diftab_offset = diftab_dim1.
    let diftab_off = diftab_dim1;
    let somtab_dim1 = *nbpnts / 2 + 1;
    // OCCT: somtab_offset = somtab_dim1.
    let somtab_off = somtab_dim1;
    // gsstab is used WITHOUT pointer adjustment (raw index ir + igss*dim1).

    let ibb = sys_base::mnfndeb_();
    if ibb >= 2 {
        sys_base::mgenmsg_("MMMAPCO");
    }
    let ikdeb = (*iordre + 1) << 1;
    let nbroot = *nbpnts / 2;

    for nd in 1..=*ndim {
        // ----------------- Calculate the coefficients of even degree ----------
        let mut ik = ikdeb;
        while ik <= *ndgjac {
            let igss = ik - ikdeb;
            let mut bidon = 0.;
            for ir in 1..=nbroot {
                // OCCT: somtab[ir + nd * somtab_dim1] - somtab_offset.
                bidon += somtab[((ir + nd * somtab_dim1) - somtab_off) as usize]
                    * gsstab[(ir + igss * gsstab_dim1) as usize];
            }
            // OCCT: crvjac[ik + nd * crvjac_dim1] - crvjac_offset.
            crvjac[((ik + nd * crvjac_dim1) - crvjac_off) as usize] = bidon;
            ik += 2;
        }

        // --------------- Calculate the coefficients of uneven degree ----------
        let mut ik = ikdeb + 1;
        while ik <= *ndgjac {
            let igss = ik - ikdeb;
            let mut bidon = 0.;
            for ir in 1..=nbroot {
                // OCCT: diftab[ir + nd * diftab_dim1] - diftab_offset.
                bidon += diftab[((ir + nd * diftab_dim1) - diftab_off) as usize]
                    * gsstab[(ir + igss * gsstab_dim1) as usize];
            }
            crvjac[((ik + nd * crvjac_dim1) - crvjac_off) as usize] = bidon;
            ik += 2;
        }
    }

    // ------- Add terms connected to the supplementary root (0.D0) ------
    if *nbpnts % 2 == 0 {
        if ibb >= 2 {
            sys_base::mgsomsg_("MMMAPCO");
        }
        return;
    }
    for nd in 1..=*ndim {
        let mut ik = ikdeb;
        while ik <= *ndgjac {
            let igss = ik - ikdeb;
            crvjac[((ik + nd * crvjac_dim1) - crvjac_off) as usize] +=
                somtab[((0 + nd * somtab_dim1) - somtab_off) as usize]
                    * gsstab[(igss * gsstab_dim1) as usize];
            ik += 2;
        }
    }

    // ------------------------------ The end -------------------------------
    if ibb >= 2 {
        sys_base::mgsomsg_("MMMAPCO");
    }
}

/// OCCT mmaperm_ (AdvApp2Var_ApproxF2var.cxx L8198-8309) - square root of the
/// average quadratic error when only the first NCFNEW coefficients are kept.
pub fn mmaperm_(
    ncofmx: &mut i32,
    ndim: &mut i32,
    ncoeff: &mut i32,
    iordre: &mut i32,
    crvjac: &[f64],
    ncfnew: &mut i32,
    errmoy: &mut f64,
) {
    let crvjac_dim1 = *ncofmx;
    // OCCT: crvjac_offset = crvjac_dim1 + 1.
    let crvjac_off = crvjac_dim1 + 1;
    let cj = |i: i32, nd: i32| -> usize {
        // OCCT: crvjac[i + nd * crvjac_dim1] - crvjac_offset.
        ((i + nd * crvjac_dim1) - crvjac_off) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 2 {
        sys_base::mgenmsg_("MMAPERM");
    }

    // --------- Minimum degree that can be reached : Stop at 1 or IA -------
    let ia = (*iordre + 1) << 1;
    let mut ncfcut = ia + 1;
    if *ncfnew + 1 > ncfcut {
        ncfcut = *ncfnew + 1;
    }

    // -------------- Elimination of coefficients of high degree ------------
    *errmoy = 0.;
    let mut bid = 0.;
    for nd in 1..=*ndim {
        for i__ in ncfcut..=*ncoeff {
            let bidj = crvjac[cj(i__, nd)];
            bid += bidj * bidj;
        }
    }

    // ----------- Square Root of average quadratic error -----------
    bid /= 2.;
    *errmoy = bid.sqrt();

    // ------------------------------- The end ------------------------------
    if ibb >= 2 {
        sys_base::mgsomsg_("MMAPERM");
    }
}

/// OCCT AdvApp2Var_ApproxF2var::mmapptt_ (AdvApp2Var_ApproxF2var.cxx
/// L8313-8641) - loads the Gauss integration coefficients for the Jacobi
/// base from the MMAPGSS/MMAPGS0/MMAPGS1/MMAPGS2 block data.  `cgauss` is
/// dimensioned (0..nbpnts/2, 0..ndgjac-2*(jordre+1)) Fortran-wise with dim1
/// = nbpnts/2 + 1; the rcad slice is the physical buffer.
pub fn mmapptt_(ndgjac: &i32, nbpnts: &i32, jordre: &i32, cgauss: &mut [f64], iercod: &mut i32) {
    let cgauss_dim1 = nbpnts / 2 + 1;

    let ibb = sys_base::mnfndeb_();
    if ibb >= 2 {
        sys_base::mgenmsg_("MMAPPTT");
    }
    *iercod = 0;

    // ------------------- Tests on the validity of inputs ----------------
    let infdg = (*jordre + 1) << 1;
    if *nbpnts != 8
        && *nbpnts != 10
        && *nbpnts != 15
        && *nbpnts != 20
        && *nbpnts != 25
        && *nbpnts != 30
        && *nbpnts != 40
        && *nbpnts != 50
        && *nbpnts != 61
    {
        *iercod = 11;
        finish(ibb, iercod);
        return;
    }

    if *jordre < -1 || *jordre > 2 {
        *iercod = 21;
        finish(ibb, iercod);
        return;
    }

    if *ndgjac >= *nbpnts || *ndgjac < infdg {
        *iercod = 31;
        finish(ibb, iercod);
        return;
    }

    // --------------- Calculation of the start pointer following NBPNTS -----------
    let mut iptdb = 0;
    if *nbpnts > 8 {
        iptdb += (8 - infdg) << 2;
    }
    if *nbpnts > 10 {
        iptdb += (10 - infdg) * 5;
    }
    if *nbpnts > 15 {
        iptdb += (15 - infdg) * 7;
    }
    if *nbpnts > 20 {
        iptdb += (20 - infdg) * 10;
    }
    if *nbpnts > 25 {
        iptdb += (25 - infdg) * 12;
    }
    if *nbpnts > 30 {
        iptdb += (30 - infdg) * 15;
    }
    if *nbpnts > 40 {
        iptdb += (40 - infdg) * 20;
    }
    if *nbpnts > 50 {
        iptdb += (50 - infdg) * 25;
    }

    let mut ipdb0 = 1;
    if *nbpnts > 15 {
        ipdb0 += (14 - infdg) / 2 + 1;
    }
    if *nbpnts > 25 {
        ipdb0 += (24 - infdg) / 2 + 1;
    }

    // ------------------ Choice of the common depending on JORDRE -------------
    match *jordre {
        -1 => {
            // ---------------- Common MMAPGSS (case without constraints) --------
            let ilong = (*nbpnts / 2) << 3;
            for kjac in 0..=*ndgjac {
                let iptt = iptdb + kjac * (*nbpnts / 2) + 1;
                // OCCT: mcrfill_(&ilong, &gslxjs[iptt - 1],
                //                &cgauss[kjac * cgauss_dim1 + 1]);
                let src = &GSLXJS[(iptt - 1) as usize..];
                let dst = (kjac * cgauss_dim1 + 1) as usize;
                sys_base::mcrfill_(ilong, src, &mut cgauss[dst..]);
            }
            // --> Case when the number of points is uneven.
            if *nbpnts % 2 == 1 {
                let mut iptt = ipdb0;
                for kjac in (0..=*ndgjac).step_by(2) {
                    // OCCT: cgauss[kjac * cgauss_dim1] = gsl0js[iptt - 1].
                    cgauss[(kjac * cgauss_dim1) as usize] = GSL0JS[(iptt - 1) as usize];
                    iptt += 1;
                }
                for kjac in (1..=*ndgjac).step_by(2) {
                    cgauss[(kjac * cgauss_dim1) as usize] = 0.;
                }
            }
        }
        0 => {
            // ---------------- Common MMAPGS0 (case with constraints C0) --------
            let mxjac = *ndgjac - infdg;
            let ilong = (*nbpnts / 2) << 3;
            for kjac in 0..=mxjac {
                let iptt = iptdb + kjac * (*nbpnts / 2) + 1;
                let src = &GSLXJ0[(iptt - 1) as usize..];
                let dst = (kjac * cgauss_dim1 + 1) as usize;
                sys_base::mcrfill_(ilong, src, &mut cgauss[dst..]);
            }
            if *nbpnts % 2 == 1 {
                let mut iptt = ipdb0;
                for kjac in (0..=mxjac).step_by(2) {
                    cgauss[(kjac * cgauss_dim1) as usize] = GSL0J0[(iptt - 1) as usize];
                    iptt += 1;
                }
                for kjac in (1..=mxjac).step_by(2) {
                    cgauss[(kjac * cgauss_dim1) as usize] = 0.;
                }
            }
        }
        1 => {
            // ---------------- Common MMAPGS1 (case with constraints C1) --------
            let mxjac = *ndgjac - infdg;
            let ilong = (*nbpnts / 2) << 3;
            for kjac in 0..=mxjac {
                let iptt = iptdb + kjac * (*nbpnts / 2) + 1;
                let src = &GSLXJ1[(iptt - 1) as usize..];
                let dst = (kjac * cgauss_dim1 + 1) as usize;
                sys_base::mcrfill_(ilong, src, &mut cgauss[dst..]);
            }
            if *nbpnts % 2 == 1 {
                let mut iptt = ipdb0;
                for kjac in (0..=mxjac).step_by(2) {
                    cgauss[(kjac * cgauss_dim1) as usize] = GSL0J1[(iptt - 1) as usize];
                    iptt += 1;
                }
                for kjac in (1..=mxjac).step_by(2) {
                    cgauss[(kjac * cgauss_dim1) as usize] = 0.;
                }
            }
        }
        _ => {
            // ---------------- Common MMAPGS2 (case with constraints C2) --------
            let mxjac = *ndgjac - infdg;
            let ilong = (*nbpnts / 2) << 3;
            for kjac in 0..=mxjac {
                let iptt = iptdb + kjac * (*nbpnts / 2) + 1;
                let src = &GSLXJ2[(iptt - 1) as usize..];
                let dst = (kjac * cgauss_dim1 + 1) as usize;
                sys_base::mcrfill_(ilong, src, &mut cgauss[dst..]);
            }
            // --> Cas of uneven number of points.
            if *nbpnts % 2 == 1 {
                let mut iptt = ipdb0;
                for kjac in (0..=mxjac).step_by(2) {
                    cgauss[(kjac * cgauss_dim1) as usize] = GSL0J2[(iptt - 1) as usize];
                    iptt += 1;
                }
                for kjac in (1..=mxjac).step_by(2) {
                    cgauss[(kjac * cgauss_dim1) as usize] = 0.;
                }
            }
        }
    }

    // -------------------------------- The end -----------------------------
    finish(ibb, iercod);
}

/// OCCT mmapptt_ L9999 tail (maermsg_ + mgsomsg_).
fn finish(ibb: i32, iercod: &mut i32) {
    if *iercod > 0 {
        sys_base::maermsg_("MMAPPTT", iercod);
    }
    if ibb >= 2 {
        sys_base::mgsomsg_("MMAPPTT");
    }
}

/// OCCT AdvApp2Var_ApproxF2var::mma2fx6_ (AdvApp2Var_ApproxF2var.cxx
/// L7285-7572) - degree reduction when the patches are constraint patches.
/// `ndimse` / `epsapr` are 1-based emulated (index 0 unused); `ncoefu` /
/// `ncoefv` / `errmax` / `patcan` / `epsfro` keep the OCCT Fortran
/// dimensioning with the virtual origin subtracted per access.
pub fn mma2fx6_(
    ncfmxu: &mut i32,
    ncfmxv: &mut i32,
    ndimen: &mut i32,
    nbsesp: &mut i32,
    ndimse: &[i32],
    nbupat: &mut i32,
    nbvpat: &mut i32,
    iordru: &mut i32,
    iordrv: &mut i32,
    epsapr: &[f64],
    epsfro: &[f64],
    patcan: &mut [f64],
    errmax: &mut [f64],
    ncoefu: &mut [i32],
    ncoefv: &mut [i32],
) {
    let epsfro_dim1 = *nbsesp;
    // OCCT: epsfro_offset = epsfro_dim1 * 5 + 1.
    let epsfro_off = epsfro_dim1 * 5 + 1;
    let ncoefu_dim1 = *nbupat;
    // OCCT: ncoefu_offset = ncoefu_dim1 + 1.
    let ncoefu_off = ncoefu_dim1 + 1;
    let ncoefv_dim1 = *nbupat;
    let ncoefv_off = ncoefv_dim1 + 1;
    let errmax_dim1 = *nbsesp;
    let errmax_dim2 = *nbupat;
    // OCCT: errmax_offset = errmax_dim1 * (errmax_dim2 + 1) + 1.
    let errmax_off = errmax_dim1 * (errmax_dim2 + 1) + 1;
    let patcan_dim1 = *ncfmxu;
    let patcan_dim2 = *ncfmxv;
    let patcan_dim3 = *ndimen;
    let patcan_dim4 = *nbupat;
    // OCCT: patcan_offset =
    //   patcan_dim1 * (patcan_dim2 * (patcan_dim3 * (patcan_dim4 + 1) + 1) + 1) + 1.
    let patcan_off =
        patcan_dim1 * (patcan_dim2 * (patcan_dim3 * (patcan_dim4 + 1) + 1) + 1) + 1;

    let nc = |tab: &[i32], ii: i32, jj: i32, dim1: i32, off: i32| -> usize {
        // OCCT: tab[ii + jj * dim1] - off.
        ((ii + jj * dim1) - off) as usize
    };
    let pc = |a: i32, b: i32, id: i32, ii: i32, jj: i32| -> usize {
        // OCCT: patcan[a + (b + (id + (ii + jj * dim4) * dim3) * dim2) * dim1] - offset.
        ((a + (b + (id + (ii + jj * patcan_dim4) * patcan_dim3) * patcan_dim2) * patcan_dim1)
            - patcan_off) as usize
    };

    let ibb = sys_base::mnfndeb_();
    if ibb >= 3 {
        sys_base::mgenmsg_("MMA2FX6");
    }

    for jj in 1..=*nbvpat {
        for ii in 1..=*nbupat {
            let mut ncfu = ncoefu[nc(ncoefu, ii, jj, ncoefu_dim1, ncoefu_off)];
            let mut ncfv = ncoefv[nc(ncoefv, ii, jj, ncoefv_dim1, ncoefv_off)];

            // -------------------- Reduction of degree by U -------------------------
            // OCCT L200 loop (goto L300 leaves the U reduction).
            'l200: loop {
                if ncfu <= (*iordru + 1) << 1 && ncfu > 2 {
                    let mut idim = 0;
                    for ns in 1..=*nbsesp {
                        let mut tol = epsapr[(ns - 1) as usize];
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 9 - epsfro_off) as usize]);
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 10 - epsfro_off) as usize]);
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 11 - epsfro_off) as usize]);
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 12 - epsfro_off) as usize]);
                        if ii == 1 || ii == *nbupat || jj == 1 || jj == *nbvpat {
                            tol = tol.min(epsfro[(ns + epsfro_dim1 * 5 - epsfro_off) as usize]);
                            tol = tol.min(epsfro[(ns + epsfro_dim1 * 6 - epsfro_off) as usize]);
                            tol = tol.min(epsfro[(ns + epsfro_dim1 * 7 - epsfro_off) as usize]);
                            tol = tol.min(epsfro[(ns + (epsfro_dim1 << 3) - epsfro_off) as usize]);
                        }
                        let mut bid = 0.;

                        for nd in 1..=ndimse[(ns - 1) as usize] {
                            let id = idim + nd;
                            for kv in 1..=ncfv {
                                // OCCT: bid += abs(patcan[ncfu + (kv + (id + ...)
                                //   * patcan_dim3) * patcan_dim2) * patcan_dim1]).
                                bid += patcan[pc(ncfu, kv, id, ii, jj)].abs();
                            }
                        }

                        if bid > tol * 1e-6
                            || bid > errmax[(ns + (ii + jj * errmax_dim2) * errmax_dim1 - errmax_off)
                                as usize]
                        {
                            break 'l200; // OCCT goto L300.
                        }
                        idim += ndimse[(ns - 1) as usize];
                    }

                    ncfu -= 1;
                    continue 'l200;
                }
                break; // OCCT fall-through to L300.
            }

            // -------------------- Reduction of degree by V -------------------------
            // OCCT L300 loop (goto L400 leaves the V reduction).
            'l300: loop {
                if ncfv <= (*iordrv + 1) << 1 && ncfv > 2 {
                    let mut idim = 0;
                    for ns in 1..=*nbsesp {
                        let mut tol = epsapr[(ns - 1) as usize];
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 9 - epsfro_off) as usize]);
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 10 - epsfro_off) as usize]);
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 11 - epsfro_off) as usize]);
                        tol = tol.min(epsfro[(ns + epsfro_dim1 * 12 - epsfro_off) as usize]);
                        if ii == 1 || ii == *nbupat || jj == 1 || jj == *nbvpat {
                            tol = tol.min(epsfro[(ns + epsfro_dim1 * 5 - epsfro_off) as usize]);
                            tol = tol.min(epsfro[(ns + epsfro_dim1 * 6 - epsfro_off) as usize]);
                            tol = tol.min(epsfro[(ns + epsfro_dim1 * 7 - epsfro_off) as usize]);
                            tol = tol.min(epsfro[(ns + (epsfro_dim1 << 3) - epsfro_off) as usize]);
                        }
                        let mut bid = 0.;

                        for nd in 1..=ndimse[(ns - 1) as usize] {
                            let id = idim + nd;
                            for ku in 1..=ncfu {
                                // OCCT: bid += abs(patcan[ku + (ncfv + (id + ...)
                                //   * patcan_dim3) * patcan_dim2) * patcan_dim1]).
                                bid += patcan[pc(ku, ncfv, id, ii, jj)].abs();
                            }
                        }

                        if bid > tol * 1e-6
                            || bid > errmax[(ns + (ii + jj * errmax_dim2) * errmax_dim1 - errmax_off)
                                as usize]
                        {
                            break 'l300; // OCCT goto L400.
                        }
                        idim += ndimse[(ns - 1) as usize];
                    }

                    ncfv -= 1;
                    continue 'l300;
                }
                break; // OCCT fall-through to L400.
            }

            // --- Return the nbs of coeff. and pass to the next square --- (L400)
            ncoefu[nc(ncoefu, ii, jj, ncoefu_dim1, ncoefu_off)] = ncfu.max(2);
            ncoefv[nc(ncoefv, ii, jj, ncoefv_dim1, ncoefv_off)] = ncfv.max(2);
        }
    }

    // ------------------------------ The End -------------------------------
    if ibb >= 3 {
        sys_base::mgsomsg_("MMA2FX6");
    }
}
