//! OCCT AdvApp2Var_Iso (AdvApp2Var_Iso.hxx + AdvApp2Var_Iso.cxx) - used to
//! store constraints on a line U = Ui or V = Vj.
//!
//! Encoding notes (architecture):
//! - `occ::handle<NCollection_HArray1/2<double>>` members ->
//!   `Option<Vec<f64>>` / `Option<Array2<f64>>` (None = null handle; the
//!   ctors leave them null, MakeApprox fills them).
//! - `GeomAbs_IsoType` -> `rcad_kernel` `IsoType` (IsoU / IsoV; the OCCT
//!   GeomAbs_NoneIso enumerator does not exist in the rcad kernel enum - the
//!   Iso ctors only accept IsoU/IsoV so the OCCT default branch is dead).
//! - the shared-handle semantics of MakeApprox (in-place mutation through
//!   the stored handle) become direct `&mut self` mutation.

use glam::DVec3;

use rcad_kernel::base::proj_lib::proj_lib_projected_curve::IsoType;

use super::approxf2var_c::EvaluatorFunc2Var;
use super::context::Context;
use super::nc_array::Array2;
use super::node::Node;

/// OCCT AdvApp2Var_Iso (AdvApp2Var_Iso.hxx L37-136).
#[derive(Debug, Clone)]
pub struct Iso {
    /// hxx L119: GeomAbs_IsoType myType.
    my_type: IsoType,
    /// hxx L120: double myConstPar.
    my_const_par: f64,
    /// hxx L121: double myU0.
    my_u0: f64,
    /// hxx L122: double myU1.
    my_u1: f64,
    /// hxx L123: double myV0.
    my_v0: f64,
    /// hxx L124: double myV1.
    my_v1: f64,
    /// hxx L125: int myPosition.
    my_position: i32,
    /// hxx L126: int myExtremOrder.
    my_extrem_order: i32,
    /// hxx L127: int myDerivOrder.
    my_deriv_order: i32,
    /// hxx L128: int myNbCoeff.
    my_nb_coeff: i32,
    /// hxx L129: bool myApprIsDone.
    my_appr_is_done: bool,
    /// hxx L130: bool myHasResult.
    my_has_result: bool,
    /// hxx L131: handle(NCollection_HArray1<double>) myEquation.
    my_equation: Option<Vec<f64>>,
    /// hxx L132: handle(NCollection_HArray2<double>) myMaxErrors.
    my_max_errors: Option<Array2<f64>>,
    /// hxx L133: handle(NCollection_HArray2<double>) myMoyErrors.
    my_moy_errors: Option<Array2<f64>>,
    /// hxx L134: handle(NCollection_HArray1<double>) mySomTab.
    my_som_tab: Option<Vec<f64>>,
    /// hxx L135: handle(NCollection_HArray1<double>) myDifTab.
    my_dif_tab: Option<Vec<f64>>,
}

impl Iso {
    /// OCCT AdvApp2Var_Iso() (AdvApp2Var_Iso.cxx L33-47).
    pub fn new() -> Self {
        Iso {
            my_type: IsoType::IsoU,
            my_const_par: 0.5,
            my_u0: 0.,
            my_u1: 1.,
            my_v0: 0.,
            my_v1: 1.,
            my_position: 0,
            my_extrem_order: 2,
            my_deriv_order: 2,
            my_nb_coeff: 0,
            my_appr_is_done: false,
            my_has_result: false,
            my_equation: None,
            my_max_errors: None,
            my_moy_errors: None,
            my_som_tab: None,
            my_dif_tab: None,
        }
    }

    /// OCCT AdvApp2Var_Iso(const GeomAbs_IsoType type, const int iu,
    /// const int iv) (AdvApp2Var_Iso.cxx L51-54) - delegates to the full
    /// ctor with (type, 0.5, 0.0, 1.0, 0.0, 1.0, 0, iu, iv).
    pub fn new_orders(the_type: IsoType, iu: i32, iv: i32) -> Self {
        Self::new_full(the_type, 0.5, 0.0, 1.0, 0.0, 1.0, 0, iu, iv)
    }

    /// OCCT AdvApp2Var_Iso(type, cte, Ufirst, Ulast, Vfirst, Vlast, pos, iu,
    /// iv) (AdvApp2Var_Iso.cxx L58-88).
    #[allow(clippy::too_many_arguments)]
    pub fn new_full(
        the_type: IsoType,
        cte: f64,
        ufirst: f64,
        ulast: f64,
        vfirst: f64,
        vlast: f64,
        pos: i32,
        iu: i32,
        iv: i32,
    ) -> Self {
        // Member init-list (L67-76); myExtremOrder/myDerivOrder are set
        // below (L78-87).
        let mut r = Iso {
            my_type: the_type,
            my_const_par: cte,
            my_u0: ufirst,
            my_u1: ulast,
            my_v0: vfirst,
            my_v1: vlast,
            my_position: pos,
            my_extrem_order: 0,
            my_deriv_order: 0,
            my_nb_coeff: 0,
            my_appr_is_done: false,
            my_has_result: false,
            my_equation: None,
            my_max_errors: None,
            my_moy_errors: None,
            my_som_tab: None,
            my_dif_tab: None,
        };
        if r.my_type == IsoType::IsoU {
            r.my_extrem_order = iv;
            r.my_deriv_order = iu;
        } else {
            r.my_extrem_order = iu;
            r.my_deriv_order = iv;
        }
        r
    }

    /// OCCT IsApproximated() (AdvApp2Var_Iso.cxx L92-95).
    pub fn is_approximated(&self) -> bool {
        self.my_appr_is_done
    }

    /// OCCT HasResult() (AdvApp2Var_Iso.cxx L99-102).
    pub fn has_result(&self) -> bool {
        self.my_has_result
    }

    /// OCCT MakeApprox(Conditions, U0, U1, V0, V1, Func, NodeBegin, NodeEnd)
    /// (AdvApp2Var_Iso.cxx L106-376).  The dead initializers of
    /// ISOFAV/NBROOT/NDGJAC/NCFLIM (= 0/0/0/1, cxx L135) are OCCT-faithful.
    #[allow(unused_assignments)]
    pub fn make_approx(
        &mut self,
        conditions: &Context,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
        func: &dyn EvaluatorFunc2Var,
        node_begin: &mut Node,
        node_end: &mut Node,
    ) {
        // fixed values
        let mut nbcrmx: i32 = 1;
        let mut nbcrbe: i32 = 0; // OCCT: uninitialized int NBCRBE (output).
        // data stored in the Context
        let mut ndimen = conditions.total_dimension();
        let mut nbsesp = conditions.total_number_ssp();
        // Attention : works only in 3D
        let ndimse: i32 = 3;
        // the domain of the grid
        let mut uvfonc = [0.0f64; 4];
        uvfonc[0] = u0;
        uvfonc[1] = u1;
        uvfonc[2] = v0;
        uvfonc[3] = v1;

        // data related to the processed iso
        let mut iordre = self.my_extrem_order;
        let ideriv = self.my_deriv_order;
        let mut tconst = self.my_const_par;

        // data related to the type of the iso
        let mut isofav: i32 = 0;
        let mut nbroot: i32 = 0;
        let mut ndgjac: i32 = 0;
        let mut ncflim: i32 = 1;
        let mut tabdec = [0.0f64; 2];
        // HUROOT / HVROOT: handle copies (clone of the Context arrays); the
        // null-handle dereference below cannot occur once the Context is
        // fully built (OCCT would crash the same way).
        let huroot = conditions.u_roots().expect("Iso::MakeApprox : URoots");
        let hvroot = conditions.v_roots().expect("Iso::MakeApprox : VRoots");
        // ROOTLG: raw base address of the selected roots array.
        let mut rootlg: &[f64] = &[];
        match self.my_type {
            IsoType::IsoV => {
                isofav = 2;
                tabdec[0] = self.my_u0;
                tabdec[1] = self.my_u1;
                uvfonc[0] = self.my_u0;
                uvfonc[1] = self.my_u1;
                nbroot = huroot.len() as i32;
                if self.my_extrem_order > -1 {
                    nbroot -= 2;
                }
                rootlg = &huroot[..];
                ndgjac = conditions.u_jac_deg();
                ncflim = conditions.u_limit();
            }
            IsoType::IsoU => {
                isofav = 1;
                tabdec[0] = self.my_v0;
                tabdec[1] = self.my_v1;
                uvfonc[2] = self.my_v0;
                uvfonc[3] = self.my_v1;
                nbroot = hvroot.len() as i32;
                if self.my_extrem_order > -1 {
                    nbroot -= 2;
                }
                rootlg = &hvroot[..];
                ndgjac = conditions.v_jac_deg();
                ncflim = conditions.v_limit();
            }
            // OCCT GeomAbs_NoneIso / default: break (dead: the rcad IsoType
            // has no NoneIso variant and the ctors only accept IsoU/IsoV).
        }

        // data relative to the position of iso (front or cut line)
        // HEPSAPR = new NCollection_HArray1<double>(1, NBSESP).
        let mut epsapr = vec![0.0f64; nbsesp as usize];
        match self.my_position {
            0 => {
                let ctol = conditions.c_toler().expect("Iso::MakeApprox : CToler");
                for iesp in 1..=nbsesp {
                    epsapr[(iesp - 1) as usize] = ctol.value(iesp, 1);
                }
            }
            1 => {
                let ftol = conditions.f_toler().expect("Iso::MakeApprox : FToler");
                for iesp in 1..=nbsesp {
                    epsapr[(iesp - 1) as usize] = ftol.value(iesp, 1);
                }
            }
            2 => {
                let ftol = conditions.f_toler().expect("Iso::MakeApprox : FToler");
                for iesp in 1..=nbsesp {
                    epsapr[(iesp - 1) as usize] = ftol.value(iesp, 2);
                }
            }
            3 => {
                let ftol = conditions.f_toler().expect("Iso::MakeApprox : FToler");
                for iesp in 1..=nbsesp {
                    epsapr[(iesp - 1) as usize] = ftol.value(iesp, 3);
                }
            }
            4 => {
                let ftol = conditions.f_toler().expect("Iso::MakeApprox : FToler");
                for iesp in 1..=nbsesp {
                    epsapr[(iesp - 1) as usize] = ftol.value(iesp, 4);
                }
            }
            _ => {}
        }

        // the tables of approximations
        let szcrb = (ndimen * ncflim) as usize;
        // HCOURBE = new NCollection_HArray1<double>(1, SZCRB * (IDERIV + 1)).
        let mut hcourbe = vec![0.0f64; szcrb * (ideriv + 1) as usize];
        // CRBAPP walks COURBE with an offset of SZCRB per derivative order.
        let mut crbapp_off: usize = 0;
        let sztab = ((1 + nbroot / 2) * ndimen) as usize;
        let mut hsomtab = vec![0.0f64; sztab * (ideriv + 1) as usize];
        let mut somapp_off: usize = 0;
        let mut hdiftab = vec![0.0f64; sztab * (ideriv + 1) as usize];
        let mut difapp_off: usize = 0;
        // HCONTR1 / HCONTR2 = new NCollection_HArray1<double>(1,
        //   (IORDRE + 2) * NDIMEN) - 1-based flat addressing.
        let mut contr1 = vec![0.0f64; ((iordre + 2) * ndimen) as usize];
        let mut contr2 = vec![0.0f64; ((iordre + 2) * ndimen) as usize];
        // HERRMAX / HERRMOY = new NCollection_HArray2<double>(1, NBSESP, 1,
        //   IDERIV + 1).
        let mut herrmax = Array2::new(1, nbsesp, 1, ideriv + 1);
        let mut emxapp = vec![0.0f64; nbsesp as usize];
        let mut herrmoy = Array2::new(1, nbsesp, 1, ideriv + 1);
        let mut emyapp = vec![0.0f64; nbsesp as usize];
        //
        // the approximations
        //
        let mut iercod: i32 = 0;
        let mut ncoeff: i32 = 0;
        // OCCT for (iapp = 0; iapp <= IDERIV; iapp++) - kept as a counted
        // loop so the engine can receive the loop variable by address
        // (OCCT passes &iapp).
        let mut iapp: i32 = 0;
        while iapp <= ideriv {
            //   approximation of the derivative of order iapp
            let mut ncfapp: i32 = 0;
            let mut ierapp: i32 = 0;
            // OCCT: mma2fnc_(&NDIMEN, &NBSESP, &NDIMSE, UVFONC, Func,
            //   &TCONST, &ISOFAV, &NBROOT, ROOTLG, &IORDRE, &iapp, &NDGJAC,
            //   &NBCRMX, &NCFLIM, EPSAPR, &ncfapp, CRBAPP, &NBCRBE, SOMAPP,
            //   DIFAPP, CONTR1, CONTR2, TABDEC, EMXAPP, EMYAPP, &ierapp).
            super::approxf2var_f::mma2fnc_(
                &mut ndimen,
                &mut nbsesp,
                std::slice::from_ref(&ndimse),
                &uvfonc,
                func,
                &mut tconst,
                &mut isofav,
                &mut nbroot,
                rootlg,
                &mut iordre,
                &mut iapp,
                &mut ndgjac,
                &mut nbcrmx,
                &mut ncflim,
                &epsapr,
                std::slice::from_mut(&mut ncfapp),
                &mut hcourbe[crbapp_off..],
                &mut nbcrbe,
                &mut hsomtab[somapp_off..],
                &mut hdiftab[difapp_off..],
                &mut contr1,
                &mut contr2,
                &mut tabdec,
                &mut emxapp,
                &mut emyapp,
                &mut ierapp,
            );
            //   error and coefficient management.
            if ierapp > 0 {
                self.my_appr_is_done = false;
                self.my_has_result = false;
                // OCCT goto FINISH (the trailing delete[] calls are no-ops
                // under Drop).
                return;
            }
            if ncoeff <= ncfapp {
                ncoeff = ncfapp;
            }
            if ierapp == -1 {
                iercod = -1;
            }
            //   return constraints of order 0 to IORDRE of extremities
            // jpos = HCONTR1->Lower() (= 1); 1-based -> vec index jpos - 1.
            let mut jpos: i32 = 1;
            for ider in 0..=iordre {
                let pt = DVec3::new(
                    contr1[(jpos - 1) as usize],
                    contr1[jpos as usize],
                    contr1[(jpos + 1) as usize],
                );
                if isofav == 2 {
                    node_begin.set_point(ider, iapp, pt);
                } else {
                    node_begin.set_point(iapp, ider, pt);
                }
                jpos += 3;
            }
            let mut jpos: i32 = 1;
            for ider in 0..=iordre {
                let pt = DVec3::new(
                    contr2[(jpos - 1) as usize],
                    contr2[jpos as usize],
                    contr2[(jpos + 1) as usize],
                );
                if isofav == 2 {
                    node_end.set_point(ider, iapp, pt);
                } else {
                    node_end.set_point(iapp, ider, pt);
                }
                jpos += 3;
            }
            //   return errors
            for iesp in 1..=nbsesp {
                herrmax.set_value(iesp, iapp + 1, emxapp[(iesp - 1) as usize]);
                herrmoy.set_value(iesp, iapp + 1, emyapp[(iesp - 1) as usize]);
            }
            // passage to the approximation of higher order
            crbapp_off += szcrb;
            somapp_off += sztab;
            difapp_off += sztab;
            iapp += 1;
        }

        // management of results
        if iercod == 0 {
            //   all approximations are correct
            self.my_appr_is_done = true;
            self.my_has_result = true;
        } else if iercod == -1 {
            //   at least one approximation is not correct
            self.my_appr_is_done = false;
            self.my_has_result = true;
        } else {
            self.my_appr_is_done = false;
            self.my_has_result = false;
        }
        if self.my_has_result {
            self.my_equation = Some(hcourbe);
            self.my_nb_coeff = ncoeff;
            self.my_max_errors = Some(herrmax);
            self.my_moy_errors = Some(herrmoy);
            self.my_som_tab = Some(hsomtab);
            self.my_dif_tab = Some(hdiftab);
        }
        // OCCT FINISH: delete[] EMXAPP / EMYAPP (no-op under Drop).
    }

    /// OCCT ChangeDomain(a, b) (AdvApp2Var_Iso.cxx L380-392).
    pub fn change_domain(&mut self, a: f64, b: f64) {
        if self.my_type == IsoType::IsoU {
            self.my_v0 = a;
            self.my_v1 = b;
        } else {
            self.my_u0 = a;
            self.my_u1 = b;
        }
    }

    /// OCCT ChangeDomain(a, b, c, d) (AdvApp2Var_Iso.cxx L396-402).
    pub fn change_domain_full(&mut self, a: f64, b: f64, c: f64, d: f64) {
        self.my_u0 = a;
        self.my_u1 = b;
        self.my_v0 = c;
        self.my_v1 = d;
    }

    /// OCCT SetConstante(newcte) (AdvApp2Var_Iso.cxx L406-409).
    pub fn set_constante(&mut self, newcte: f64) {
        self.my_const_par = newcte;
    }

    /// OCCT SetPosition(newpos) (AdvApp2Var_Iso.cxx L413-416).
    pub fn set_position(&mut self, newpos: i32) {
        self.my_position = newpos;
    }

    /// OCCT ResetApprox() (AdvApp2Var_Iso.cxx L420-424).
    pub fn reset_approx(&mut self) {
        self.my_appr_is_done = false;
        self.my_has_result = false;
    }

    /// OCCT OverwriteApprox() (AdvApp2Var_Iso.cxx L428-434).
    pub fn overwrite_approx(&mut self) {
        if self.my_has_result {
            self.my_appr_is_done = true;
        }
    }

    /// OCCT Type() (AdvApp2Var_Iso.cxx L438-441).
    pub fn type_(&self) -> IsoType {
        self.my_type
    }

    /// OCCT Constante() (AdvApp2Var_Iso.cxx L445-448).
    pub fn constante(&self) -> f64 {
        self.my_const_par
    }

    /// OCCT T0() (AdvApp2Var_Iso.cxx L452-462).
    pub fn t0(&self) -> f64 {
        if self.my_type == IsoType::IsoU {
            self.my_v0
        } else {
            self.my_u0
        }
    }

    /// OCCT T1() (AdvApp2Var_Iso.cxx L466-476).
    pub fn t1(&self) -> f64 {
        if self.my_type == IsoType::IsoU {
            self.my_v1
        } else {
            self.my_u1
        }
    }

    /// OCCT U0() (AdvApp2Var_Iso.cxx L480-483).
    pub fn u0(&self) -> f64 {
        self.my_u0
    }

    /// OCCT U1() (AdvApp2Var_Iso.cxx L487-490).
    pub fn u1(&self) -> f64 {
        self.my_u1
    }

    /// OCCT V0() (AdvApp2Var_Iso.cxx L494-501).
    pub fn v0(&self) -> f64 {
        self.my_v0
    }

    /// OCCT V1() (AdvApp2Var_Iso.cxx L505-508).
    pub fn v1(&self) -> f64 {
        self.my_v1
    }

    /// OCCT UOrder() (AdvApp2Var_Iso.cxx L508-518).
    pub fn u_order(&self) -> i32 {
        if self.type_() == IsoType::IsoU {
            self.my_deriv_order
        } else {
            self.my_extrem_order
        }
    }

    /// OCCT VOrder() (AdvApp2Var_Iso.cxx L522-532).
    pub fn v_order(&self) -> i32 {
        if self.type_() == IsoType::IsoV {
            self.my_deriv_order
        } else {
            self.my_extrem_order
        }
    }

    /// OCCT Position() (AdvApp2Var_Iso.cxx L536-539).
    pub fn position(&self) -> i32 {
        self.my_position
    }

    /// OCCT NbCoeff() (AdvApp2Var_Iso.cxx L543-546).
    pub fn nb_coeff(&self) -> i32 {
        self.my_nb_coeff
    }

    /// OCCT Polynom() (AdvApp2Var_Iso.cxx L550-553) - const handle reference
    /// returned by value-copy (clone).
    pub fn polynom(&self) -> Option<Vec<f64>> {
        self.my_equation.clone()
    }

    /// OCCT SomTab() (AdvApp2Var_Iso.cxx L555-558).
    pub fn som_tab(&self) -> Option<Vec<f64>> {
        self.my_som_tab.clone()
    }

    /// OCCT DifTab() (AdvApp2Var_Iso.cxx L560-563).
    pub fn dif_tab(&self) -> Option<Vec<f64>> {
        self.my_dif_tab.clone()
    }

    /// OCCT MaxErrors() (AdvApp2Var_Iso.cxx L565-568).
    pub fn max_errors(&self) -> Option<Array2<f64>> {
        self.my_max_errors.clone()
    }

    /// OCCT MoyErrors() (AdvApp2Var_Iso.cxx L570-573).
    pub fn moy_errors(&self) -> Option<Array2<f64>> {
        self.my_moy_errors.clone()
    }
}
