//! OCCT AdvApp2Var_Context (AdvApp2Var_Context.hxx + AdvApp2Var_Context.cxx)
//! - contains all the parameters for approximation (tolerancy, computing
//! option, ...).
//!
//! Encoding notes (architecture):
//! - `occ::handle<NCollection_HArray1<double>>` members ->
//!   `Option<Vec<f64>>` (None = null handle; the default ctor leaves every
//!   handle null, the full ctor fills them).
//! - `occ::handle<NCollection_HArray2<double>>` members ->
//!   `Option<super::nc_array::Array2<f64>>`.
//! - `throw Standard_ConstructionError(...)` -> `panic!(...)`.

use super::approxf2var_b::{mma2jmx_, mma2roo_, mmapptt_};
use super::nc_array::Array2;

/// OCCT lesparam (AdvApp2Var_Context.cxx L22-78) - calculation of
/// parameters.  Returns false on the OCCT failure path (icodeo < 0).
fn lesparam(iordre: i32, ncflim: i32, icodeo: i32, nbpnts: &mut i32, ndgjac: &mut i32) -> bool {
    // jacobi degree
    *ndgjac = ncflim; // it always keeps a reserve coefficient
    if icodeo < 0 {
        return false;
    }
    if icodeo > 0 {
        *ndgjac += 9 - (iordre + 1); // iordre rescales the frequences upwards
        *ndgjac += (icodeo - 1) * 10;
    }
    // ---> Min Number of required points.
    if *ndgjac < 8 {
        *nbpnts = 8;
    } else if *ndgjac < 10 {
        *nbpnts = 10;
    }
    //  else if (ndgjac < 15) { nbpnt = 15; } Bug Uneven number
    else if *ndgjac < 20 {
        *nbpnts = 20;
    }
    //  else if (ndgjac < 25) { nbpnt = 25; } Bug Uneven number
    else if *ndgjac < 30 {
        *nbpnts = 30;
    } else if *ndgjac < 40 {
        *nbpnts = 40;
    } else if *ndgjac < 50 {
        *nbpnts = 50;
    }
    //  else if (*ndgjac < 61) { nbpnt = 61;} Bug Uneven number
    else {
        *nbpnts = 50;
        // OCCT_DEBUG print dropped (compiled out).
    }

    // If constraints are on borders, this adds 2 points
    if iordre > -1 {
        *nbpnts += 2;
    }

    true
}

/// OCCT appendInternalTolerance (AdvApp2Var_Context.cxx L80-89).  The
/// HArray1 carrier is a flat Vec addressed 1-based (index k -> [k - 1]).
fn append_internal_tolerance(
    the_source: &[f64],
    the_count: i32,
    the_offset: i32,
    the_target: &mut Vec<f64>,
) {
    for an_ssp_index in 1..=the_count {
        the_target[(the_offset + an_ssp_index - 1) as usize] =
            the_source[(an_ssp_index - 1) as usize];
    }
}

/// OCCT appendFrontierTolerance (AdvApp2Var_Context.cxx L91-106).
fn append_frontier_tolerance(
    the_source: &Array2<f64>,
    the_count: i32,
    the_offset: i32,
    the_frontier_tol: &mut Array2<f64>,
    the_cutting_tol: &mut Array2<f64>,
) {
    for an_ssp_index in 1..=the_count {
        let a_global_ssp_index = the_offset + an_ssp_index;
        for a_tol_index in 1..=4 {
            the_frontier_tol.set_value(
                a_global_ssp_index,
                a_tol_index,
                the_source.value(an_ssp_index, a_tol_index),
            );
            the_cutting_tol.set_value(a_global_ssp_index, a_tol_index, 0.0);
        }
    }
}

/// OCCT hMaxFactor (AdvApp2Var_Context.cxx L108-121).
fn h_max_factor(the_order: i32) -> f64 {
    match the_order {
        -1 => 0.0,
        0 => 1.0,
        1 => 1.5,
        _ => 1.75,
    }
}

/// OCCT AdvApp2Var_Context (AdvApp2Var_Context.hxx L32-113).
#[derive(Debug, Clone)]
pub struct Context {
    /// hxx L92: int myFav.
    my_fav: i32,
    /// hxx L93: int myOrdU.
    my_ord_u: i32,
    /// hxx L94: int myOrdV.
    my_ord_v: i32,
    /// hxx L95: int myLimU.
    my_lim_u: i32,
    /// hxx L96: int myLimV.
    my_lim_v: i32,
    /// hxx L97: int myNb1DSS.
    my_nb1_dss: i32,
    /// hxx L98: int myNb2DSS.
    my_nb2_dss: i32,
    /// hxx L99: int myNb3DSS.
    my_nb3_dss: i32,
    /// hxx L100: int myNbURoot.
    my_nb_u_root: i32,
    /// hxx L101: int myNbVRoot.
    my_nb_v_root: i32,
    /// hxx L102: int myJDegU.
    my_j_deg_u: i32,
    /// hxx L103: int myJDegV.
    my_j_deg_v: i32,
    /// hxx L104: handle(NCollection_HArray1<double>) myJMaxU.
    my_j_max_u: Option<Vec<f64>>,
    /// hxx L105: handle(NCollection_HArray1<double>) myJMaxV.
    my_j_max_v: Option<Vec<f64>>,
    /// hxx L106: handle(NCollection_HArray1<double>) myURoots.
    my_u_roots: Option<Vec<f64>>,
    /// hxx L107: handle(NCollection_HArray1<double>) myVRoots.
    my_v_roots: Option<Vec<f64>>,
    /// hxx L108: handle(NCollection_HArray1<double>) myUGauss.
    my_u_gauss: Option<Vec<f64>>,
    /// hxx L109: handle(NCollection_HArray1<double>) myVGauss.
    my_v_gauss: Option<Vec<f64>>,
    /// hxx L110: handle(NCollection_HArray1<double>) myInternalTol.
    my_internal_tol: Option<Vec<f64>>,
    /// hxx L111: handle(NCollection_HArray2<double>) myFrontierTol.
    my_frontier_tol: Option<Array2<f64>>,
    /// hxx L112: handle(NCollection_HArray2<double>) myCuttingTol.
    my_cutting_tol: Option<Array2<f64>>,
}

impl Context {
    /// OCCT AdvApp2Var_Context() (AdvApp2Var_Context.cxx L125-139) - every
    /// scalar zero, every handle null.
    pub fn new() -> Self {
        Context {
            my_fav: 0,
            my_ord_u: 0,
            my_ord_v: 0,
            my_lim_u: 0,
            my_lim_v: 0,
            my_nb1_dss: 0,
            my_nb2_dss: 0,
            my_nb3_dss: 0,
            my_nb_u_root: 0,
            my_nb_v_root: 0,
            my_j_deg_u: 0,
            my_j_deg_v: 0,
            my_j_max_u: None,
            my_j_max_v: None,
            my_u_roots: None,
            my_v_roots: None,
            my_u_gauss: None,
            my_v_gauss: None,
            my_internal_tol: None,
            my_frontier_tol: None,
            my_cutting_tol: None,
        }
    }

    /// OCCT AdvApp2Var_Context(ifav, iu, iv, nlimu, nlimv, iprecis, nb1Dss,
    /// nb2Dss, nb3Dss, tol1D, tol2D, tol3D, tof1D, tof2D, tof3D)
    /// (AdvApp2Var_Context.cxx L143-291).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_params(
        ifav: i32,
        iu: i32,
        iv: i32,
        nlimu: i32,
        nlimv: i32,
        iprecis: i32,
        nb1_dss: i32,
        nb2_dss: i32,
        nb3_dss: i32,
        tol1_d: &[f64],
        tol2_d: &[f64],
        tol3_d: &[f64],
        tof1_d: &Array2<f64>,
        tof2_d: &Array2<f64>,
        tof3_d: &Array2<f64>,
    ) -> Self {
        // Member init-list (L158-165); the handles stay null until filled.
        let mut r = Context {
            my_fav: ifav,
            my_ord_u: iu,
            my_ord_v: iv,
            my_lim_u: nlimu,
            my_lim_v: nlimv,
            my_nb1_dss: nb1_dss,
            my_nb2_dss: nb2_dss,
            my_nb3_dss: nb3_dss,
            my_nb_u_root: 0,
            my_nb_v_root: 0,
            my_j_deg_u: 0,
            my_j_deg_v: 0,
            my_j_max_u: None,
            my_j_max_v: None,
            my_u_roots: None,
            my_v_roots: None,
            my_u_gauss: None,
            my_v_gauss: None,
            my_internal_tol: None,
            my_frontier_tol: None,
            my_cutting_tol: None,
        };

        let mut an_error_code: i32 = 0;
        let mut nb_pnt_u: i32 = 0;
        let mut j_deg_u: i32 = 0;
        let mut nb_pnt_v: i32 = 0;
        let mut j_deg_v: i32 = 0;
        let mut an_order_u = iu;
        let mut an_order_v = iv;

        // myNbURoot,myJDegU
        let mut a_coeff_limit = nlimu;
        if a_coeff_limit < 2 * iu + 2 {
            a_coeff_limit = 2 * iu + 2;
        }
        if !lesparam(iu, a_coeff_limit, iprecis, &mut nb_pnt_u, &mut j_deg_u) {
            panic!("AdvApp2Var_Context");
        }
        r.my_nb_u_root = nb_pnt_u;
        r.my_j_deg_u = j_deg_u;
        if iu > -1 {
            nb_pnt_u = r.my_nb_u_root - 2;
        }

        // myJMaxU
        let a_size = (j_deg_u - 2 * iu - 1) as usize;
        // new NCollection_HArray1<double>(1, aSize); JU_array = raw base.
        let mut j_max_u = vec![0.0f64; a_size];
        // AdvApp2Var_ApproxF2var::mma2jmx_(&JDegU, &anOrderU, JU_array);
        mma2jmx_(&mut j_deg_u, &mut an_order_u, &mut j_max_u);
        r.my_j_max_u = Some(j_max_u);

        // myNbVRoot,myJDegV
        let mut a_coeff_limit = nlimv;
        if a_coeff_limit < 2 * iv + 2 {
            a_coeff_limit = 2 * iv + 2;
        }
        if !lesparam(iv, a_coeff_limit, iprecis, &mut nb_pnt_v, &mut j_deg_v) {
            panic!("AdvApp2Var_Context");
        }
        r.my_nb_v_root = nb_pnt_v;
        r.my_j_deg_v = j_deg_v;
        if iv > -1 {
            nb_pnt_v = r.my_nb_v_root - 2;
        }

        // myJMaxV
        let a_size = (j_deg_v - 2 * iv - 1) as usize;
        let mut j_max_v = vec![0.0f64; a_size];
        // AdvApp2Var_ApproxF2var::mma2jmx_(&JDegV, &anOrderV, JV_array);
        mma2jmx_(&mut j_deg_v, &mut an_order_v, &mut j_max_v);
        r.my_j_max_v = Some(j_max_v);

        // myURoots, myVRoots
        let mut u_roots = vec![0.0f64; r.my_nb_u_root as usize];
        let mut v_roots = vec![0.0f64; r.my_nb_v_root as usize];
        // AdvApp2Var_ApproxF2var::mma2roo_(&NbPntU, &NbPntV, U_array, V_array);
        mma2roo_(&mut nb_pnt_u, &mut nb_pnt_v, &mut u_roots, &mut v_roots);
        r.my_u_roots = Some(u_roots);
        r.my_v_roots = Some(v_roots);

        // myUGauss
        let a_size = ((nb_pnt_u / 2 + 1) * (r.my_j_deg_u - 2 * iu - 1)) as usize;
        let mut u_gauss = vec![0.0f64; a_size];
        // AdvApp2Var_ApproxF2var::mmapptt_(&JDegU, &NbPntU, &anOrderU,
        //                                  UG_array, &anErrorCode);
        mmapptt_(&j_deg_u, &nb_pnt_u, &an_order_u, &mut u_gauss, &mut an_error_code);
        if an_error_code != 0 {
            panic!("AdvApp2Var_Context : Error in FORTRAN");
        }
        r.my_u_gauss = Some(u_gauss);

        // myVGauss
        let a_size = ((nb_pnt_v / 2 + 1) * (r.my_j_deg_v - 2 * iv - 1)) as usize;
        let mut v_gauss = vec![0.0f64; a_size];
        // AdvApp2Var_ApproxF2var::mmapptt_(&JDegV, &NbPntV, &anOrderV,
        //                                  VG_array, &anErrorCode);
        mmapptt_(&j_deg_v, &nb_pnt_v, &an_order_v, &mut v_gauss, &mut an_error_code);
        if an_error_code != 0 {
            panic!("AdvApp2Var_Context : Error in FORTRAN");
        }
        r.my_v_gauss = Some(v_gauss);

        // myInternalTol, myFrontierTol, myCuttingTol
        let a_nb_ssp = nb1_dss + nb2_dss + nb3_dss;
        // new NCollection_HArray1<double>(1, aNbSSP) - 1-based flat Vec.
        let mut i_tol = vec![0.0f64; a_nb_ssp as usize];
        append_internal_tolerance(tol1_d, nb1_dss, 0, &mut i_tol);
        append_internal_tolerance(tol2_d, nb2_dss, nb1_dss, &mut i_tol);
        append_internal_tolerance(tol3_d, nb3_dss, nb1_dss + nb2_dss, &mut i_tol);
        if iu > -1 || iv > -1 {
            for an_ssp_index in 1..=a_nb_ssp {
                i_tol[(an_ssp_index - 1) as usize] /= 2.0;
            }
        }
        let mut f_tol = Array2::new(1, a_nb_ssp, 1, 4);
        let mut c_tol = Array2::new(1, a_nb_ssp, 1, 4);
        append_frontier_tolerance(tof1_d, nb1_dss, 0, &mut f_tol, &mut c_tol);
        append_frontier_tolerance(tof2_d, nb2_dss, nb1_dss, &mut f_tol, &mut c_tol);
        append_frontier_tolerance(tof3_d, nb3_dss, nb1_dss + nb2_dss, &mut f_tol, &mut c_tol);
        if iu > -1 || iv > -1 {
            let h_max_u = h_max_factor(iu);
            let h_max_v = h_max_factor(iv);
            let a_weight = h_max_u * h_max_v + h_max_u + h_max_v;
            for an_ssp_index in 1..=a_nb_ssp {
                for a_tol_index in 1..=4 {
                    let a_tol_min = i_tol[(an_ssp_index - 1) as usize] / a_weight;
                    if a_tol_min < f_tol.value(an_ssp_index, a_tol_index) {
                        f_tol.set_value(an_ssp_index, a_tol_index, a_tol_min);
                    }
                    c_tol.set_value(an_ssp_index, a_tol_index, a_tol_min);
                }
            }
        }
        r.my_internal_tol = Some(i_tol);
        r.my_frontier_tol = Some(f_tol);
        r.my_cutting_tol = Some(c_tol);
        r
    }

    /// OCCT TotalDimension() (AdvApp2Var_Context.cxx L295-298).
    pub fn total_dimension(&self) -> i32 {
        self.my_nb1_dss + 2 * self.my_nb2_dss + 3 * self.my_nb3_dss
    }

    /// OCCT TotalNumberSSP() (AdvApp2Var_Context.cxx L302-305).
    pub fn total_number_ssp(&self) -> i32 {
        self.my_nb1_dss + self.my_nb2_dss + self.my_nb3_dss
    }

    /// OCCT FavorIso() (AdvApp2Var_Context.cxx L309-312).
    pub fn favor_iso(&self) -> i32 {
        self.my_fav
    }

    /// OCCT UOrder() (AdvApp2Var_Context.cxx L316-319).
    pub fn u_order(&self) -> i32 {
        self.my_ord_u
    }

    /// OCCT VOrder() (AdvApp2Var_Context.cxx L323-326).
    pub fn v_order(&self) -> i32 {
        self.my_ord_v
    }

    /// OCCT ULimit() (AdvApp2Var_Context.cxx L330-333).
    pub fn u_limit(&self) -> i32 {
        self.my_lim_u
    }

    /// OCCT VLimit() (AdvApp2Var_Context.cxx L336-340).
    pub fn v_limit(&self) -> i32 {
        self.my_lim_v
    }

    /// OCCT UJacDeg() (AdvApp2Var_Context.cxx L344-347).
    pub fn u_jac_deg(&self) -> i32 {
        self.my_j_deg_u
    }

    /// OCCT VJacDeg() (AdvApp2Var_Context.cxx L351-354).
    pub fn v_jac_deg(&self) -> i32 {
        self.my_j_deg_v
    }

    /// OCCT UJacMax() (AdvApp2Var_Context.cxx L358-361) - handle copy
    /// becomes a clone of the array.
    pub fn u_jac_max(&self) -> Option<Vec<f64>> {
        self.my_j_max_u.clone()
    }

    /// OCCT VJacMax() (AdvApp2Var_Context.cxx L365-368).
    pub fn v_jac_max(&self) -> Option<Vec<f64>> {
        self.my_j_max_v.clone()
    }

    /// OCCT URoots() (AdvApp2Var_Context.cxx L372-375).
    pub fn u_roots(&self) -> Option<Vec<f64>> {
        self.my_u_roots.clone()
    }

    /// OCCT VRoots() (AdvApp2Var_Context.cxx L379-382).
    pub fn v_roots(&self) -> Option<Vec<f64>> {
        self.my_v_roots.clone()
    }

    /// OCCT UGauss() (AdvApp2Var_Context.cxx L386-389).
    pub fn u_gauss(&self) -> Option<Vec<f64>> {
        self.my_u_gauss.clone()
    }

    /// OCCT VGauss() (AdvApp2Var_Context.cxx L393-396).
    pub fn v_gauss(&self) -> Option<Vec<f64>> {
        self.my_v_gauss.clone()
    }

    /// OCCT IToler() (AdvApp2Var_Context.cxx L400-403).
    pub fn i_toler(&self) -> Option<Vec<f64>> {
        self.my_internal_tol.clone()
    }

    /// OCCT FToler() (AdvApp2Var_Context.cxx L407-410).
    pub fn f_toler(&self) -> Option<Array2<f64>> {
        self.my_frontier_tol.clone()
    }

    /// OCCT CToler() (AdvApp2Var_Context.cxx L414-417).
    pub fn c_toler(&self) -> Option<Array2<f64>> {
        self.my_cutting_tol.clone()
    }
}
