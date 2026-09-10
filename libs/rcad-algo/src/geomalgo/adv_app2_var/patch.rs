//! OCCT AdvApp2Var_Patch (AdvApp2Var_Patch.hxx + AdvApp2Var_Patch.cxx) -
//! used to store results on a domain [Ui,Ui+1]x[Vj,Vj+1].
//!
//! Encoding notes (architecture):
//! - `occ::handle<NCollection_HArray1<double>>` members ->
//!   `Option<Vec<f64>>`; `handle<NCollection_HArray2<double>>` ->
//!   `Option<Array2<f64>>` (None = null handle).
//! - `throw Standard_ConstructionError(...)` -> `panic!(...)`.
//! - Fortran engine pointers that walk a buffer (`PATCAN`, `CRBAPP`, `Is`)
//!   become base + offset sub-slices of the same Vec.

use glam::DVec3;

use super::context::Context;
use super::criterion::{Criterion, CriterionType};
use super::framework::Framework;
use super::nc_array::Array2;

/// OCCT AdvApp2Var_Patch (AdvApp2Var_Patch.hxx L37-141).
#[derive(Debug, Clone)]
pub struct Patch {
    /// hxx L120: double myU0.
    my_u0: f64,
    /// hxx L121: double myU1.
    my_u1: f64,
    /// hxx L122: double myV0.
    my_v0: f64,
    /// hxx L123: double myV1.
    my_v1: f64,
    /// hxx L124: int myOrdInU.
    my_ord_in_u: i32,
    /// hxx L125: int myOrdInV.
    my_ord_in_v: i32,
    /// hxx L126: int myNbCoeffInU.
    my_nb_coeff_in_u: i32,
    /// hxx L127: int myNbCoeffInV.
    my_nb_coeff_in_v: i32,
    /// hxx L128: bool myApprIsDone.
    my_appr_is_done: bool,
    /// hxx L129: bool myHasResult.
    my_has_result: bool,
    /// hxx L130: handle(NCollection_HArray1<double>) myEquation.
    my_equation: Option<Vec<f64>>,
    /// hxx L131: handle(NCollection_HArray1<double>) myMaxErrors.
    my_max_errors: Option<Vec<f64>>,
    /// hxx L132: handle(NCollection_HArray1<double>) myMoyErrors.
    my_moy_errors: Option<Vec<f64>>,
    /// hxx L133: handle(NCollection_HArray2<double>) myIsoErrors.
    my_iso_errors: Option<Array2<f64>>,
    /// hxx L134: int myCutSense.
    my_cut_sense: i32,
    /// hxx L135: bool myDiscIsDone.
    my_disc_is_done: bool,
    /// hxx L136: handle(NCollection_HArray1<double>) mySosoTab.
    my_soso_tab: Option<Vec<f64>>,
    /// hxx L137: handle(NCollection_HArray1<double>) myDisoTab.
    my_diso_tab: Option<Vec<f64>>,
    /// hxx L138: handle(NCollection_HArray1<double>) mySodiTab.
    my_sodi_tab: Option<Vec<f64>>,
    /// hxx L139: handle(NCollection_HArray1<double>) myDidiTab.
    my_didi_tab: Option<Vec<f64>>,
    /// hxx L140: double myCritValue.
    my_crit_value: f64,
}

impl Patch {
    /// OCCT AdvApp2Var_Patch() (AdvApp2Var_Patch.cxx L46-61).
    pub fn new() -> Self {
        Patch {
            my_u0: 0.,
            my_u1: 1.,
            my_v0: 0.,
            my_v1: 1.,
            my_ord_in_u: 0,
            my_ord_in_v: 0,
            my_nb_coeff_in_u: 0,
            my_nb_coeff_in_v: 0,
            my_appr_is_done: false,
            my_has_result: false,
            my_equation: None,
            my_max_errors: None,
            my_moy_errors: None,
            my_iso_errors: None,
            my_cut_sense: 0,
            my_disc_is_done: false,
            my_soso_tab: None,
            my_diso_tab: None,
            my_sodi_tab: None,
            my_didi_tab: None,
            my_crit_value: 0.,
        }
    }

    /// OCCT AdvApp2Var_Patch(U0, U1, V0, V1, iu, iv)
    /// (AdvApp2Var_Patch.cxx L65-85).
    pub fn new_with_domain(u0: f64, u1: f64, v0: f64, v1: f64, iu: i32, iv: i32) -> Self {
        Patch {
            my_u0: u0,
            my_u1: u1,
            my_v0: v0,
            my_v1: v1,
            my_ord_in_u: iu,
            my_ord_in_v: iv,
            my_nb_coeff_in_u: 0,
            my_nb_coeff_in_v: 0,
            my_appr_is_done: false,
            my_has_result: false,
            my_equation: None,
            my_max_errors: None,
            my_moy_errors: None,
            my_iso_errors: None,
            my_cut_sense: 0,
            my_disc_is_done: false,
            my_soso_tab: None,
            my_diso_tab: None,
            my_sodi_tab: None,
            my_didi_tab: None,
            my_crit_value: 0.,
        }
    }

    /// OCCT IsDiscretised() (AdvApp2Var_Patch.cxx L89-92).
    pub fn is_discretised(&self) -> bool {
        self.my_disc_is_done
    }

    /// OCCT Discretise(Conditions, Constraints, Func)
    /// (AdvApp2Var_Patch.cxx L96-370).
    pub fn discretise(
        &mut self,
        conditions: &Context,
        constraints: &Framework,
        func: &dyn super::approxf2var_c::EvaluatorFunc2Var,
    ) {
        // data stored in the Context
        let mut ndimen = conditions.total_dimension();
        // Attention : works only for 3D
        let mut isofav = conditions.favor_iso();

        // data related to the patch to be discretized
        let mut nbpntu: i32;
        let mut nbpntv: i32;
        let mut iordru = self.my_ord_in_u;
        let mut iordrv = self.my_ord_in_v;
        // HUROOT / HVROOT handle copies; UROOT / VROOT are the raw bases.
        let huroot = conditions.u_roots().expect("Patch::Discretise : URoots");
        let hvroot = conditions.v_roots().expect("Patch::Discretise : VRoots");
        let uroot: &[f64] = &huroot;
        nbpntu = huroot.len() as i32;
        if self.my_ord_in_u > -1 {
            nbpntu -= 2;
        }
        let vroot: &[f64] = &hvroot;
        nbpntv = hvroot.len() as i32;
        if self.my_ord_in_v > -1 {
            nbpntv -= 2;
        }

        // data stored in the Framework Constraints cad Nodes and Isos
        // C1, C2, C3 and C4 are dimensionnes in FORTRAN with
        // (NDIMEN,IORDRU+2,IORDRV+2)
        let mut size = (ndimen * (iordru + 2) * (iordrv + 2)) as usize;
        // HCOINS = new NCollection_HArray1<double>(1, SIZE * 4); Init(0.).
        let mut hcoins = vec![0.0f64; size * 4];

        let du = (self.my_u1 - self.my_u0) / 2.0;
        let dv = (self.my_v1 - self.my_v0) / 2.0;
        let mut rho: f64;
        let mut valnorm: f64;

        let mut iu: i32 = 0;
        while iu <= self.my_ord_in_u {
            let mut iv: i32 = 0;
            while iv <= self.my_ord_in_v {
                // factor of normalization
                rho = du.powi(iu) * dv.powi(iv);

                // F(U0,V0) and its derivatives normalized on (-1,1)
                valnorm = rho * constraints.node(self.my_u0, self.my_v0).point(iu, iv).x;
                hcoins[(1 + ndimen * iu + ndimen * (iordru + 2) * iv - 1) as usize] = valnorm;
                valnorm = rho * constraints.node(self.my_u0, self.my_v0).point(iu, iv).y;
                hcoins[(2 + ndimen * iu + ndimen * (iordru + 2) * iv - 1) as usize] = valnorm;
                valnorm = rho * constraints.node(self.my_u0, self.my_v0).point(iu, iv).z;
                hcoins[(3 + ndimen * iu + ndimen * (iordru + 2) * iv - 1) as usize] = valnorm;

                // F(U1,V0) and its derivatives normalized on (-1,1)
                valnorm = rho * constraints.node(self.my_u1, self.my_v0).point(iu, iv).x;
                hcoins[(size as i32 + 1 + ndimen * iu + ndimen * (iordru + 2) * iv - 1) as usize] =
                    valnorm;
                valnorm = rho * constraints.node(self.my_u1, self.my_v0).point(iu, iv).y;
                hcoins[(size as i32 + 2 + ndimen * iu + ndimen * (iordru + 2) * iv - 1) as usize] =
                    valnorm;
                valnorm = rho * constraints.node(self.my_u1, self.my_v0).point(iu, iv).z;
                hcoins[(size as i32 + 3 + ndimen * iu + ndimen * (iordru + 2) * iv - 1) as usize] =
                    valnorm;

                // F(U0,V1) and its derivatives normalized on (-1,1)
                valnorm = rho * constraints.node(self.my_u0, self.my_v1).point(iu, iv).x;
                hcoins[(2 * size as i32 + 1 + ndimen * iu + ndimen * (iordru + 2) * iv - 1)
                    as usize] = valnorm;
                valnorm = rho * constraints.node(self.my_u0, self.my_v1).point(iu, iv).y;
                hcoins[(2 * size as i32 + 2 + ndimen * iu + ndimen * (iordru + 2) * iv - 1)
                    as usize] = valnorm;
                valnorm = rho * constraints.node(self.my_u0, self.my_v1).point(iu, iv).z;
                hcoins[(2 * size as i32 + 3 + ndimen * iu + ndimen * (iordru + 2) * iv - 1)
                    as usize] = valnorm;

                // F(U1,V1) and its derivatives normalized on (-1,1)
                valnorm = rho * constraints.node(self.my_u1, self.my_v1).point(iu, iv).x;
                hcoins[(3 * size as i32 + 1 + ndimen * iu + ndimen * (iordru + 2) * iv - 1)
                    as usize] = valnorm;
                valnorm = rho * constraints.node(self.my_u1, self.my_v1).point(iu, iv).y;
                hcoins[(3 * size as i32 + 2 + ndimen * iu + ndimen * (iordru + 2) * iv - 1)
                    as usize] = valnorm;
                valnorm = rho * constraints.node(self.my_u1, self.my_v1).point(iu, iv).z;
                hcoins[(3 * size as i32 + 3 + ndimen * iu + ndimen * (iordru + 2) * iv - 1)
                    as usize] = valnorm;
                iv += 1;
            }
            iu += 1;
        }
        // C1 = base; C2 = C1 + SIZE; C3 = C2 + SIZE; C4 = C3 + SIZE.
        let c1: &[f64] = &hcoins[0..size];
        let c2: &[f64] = &hcoins[size..2 * size];
        let c3: &[f64] = &hcoins[2 * size..3 * size];
        let c4: &[f64] = &hcoins[3 * size..4 * size];

        // tables SomTab and Diftab of discretization of isos U=U0 and U=U1
        // SU0, SU1, DU0 and DU1 are dimensioned in FORTRAN to
        // (1+NBPNTV/2)*NDIMEN*(IORDRU+1)

        size = ((1 + nbpntv / 2) * ndimen) as usize;

        // HSU0 = ...; HSU0->ChangeArray1() =
        //   ((Constraints.IsoU(myU0,myV0,myV1)).SomTab())->Array1().
        let iso_u0 = constraints.iso_u(self.my_u0, self.my_v0, self.my_v1);
        let mut su0 = vec![0.0f64; size * (iordru + 1) as usize];
        su0.copy_from_slice(&iso_u0.som_tab().expect("Patch::Discretise : SomTab"));
        let iso_u0 = constraints.iso_u(self.my_u0, self.my_v0, self.my_v1);
        let mut du0 = vec![0.0f64; size * (iordru + 1) as usize];
        du0.copy_from_slice(&iso_u0.dif_tab().expect("Patch::Discretise : DifTab"));
        let iso_u1 = constraints.iso_u(self.my_u1, self.my_v0, self.my_v1);
        let mut su1 = vec![0.0f64; size * (iordru + 1) as usize];
        su1.copy_from_slice(&iso_u1.som_tab().expect("Patch::Discretise : SomTab"));
        let iso_u1 = constraints.iso_u(self.my_u1, self.my_v0, self.my_v1);
        let mut du1 = vec![0.0f64; size * (iordru + 1) as usize];
        du1.copy_from_slice(&iso_u1.dif_tab().expect("Patch::Discretise : DifTab"));

        // normalization
        let mut ideb1: i32;
        let mut ideb2: i32;
        let mut ideb3: i32;
        let mut ideb4: i32;
        let mut jj: i32;
        iu = 1;
        while iu <= iordru {
            rho = du.powi(iu);
            // OCCT: idebX = Lower + iu * SIZE - 1 (= iu * SIZE for Lower=1);
            // entries idebX + jj (1-based) -> [idebX + jj - 1].
            ideb1 = iu * size as i32;
            ideb2 = iu * size as i32;
            ideb3 = iu * size as i32;
            ideb4 = iu * size as i32;
            jj = 1;
            while jj <= size as i32 {
                su0[(ideb1 + jj - 1) as usize] = rho * su0[(ideb1 + jj - 1) as usize];
                du0[(ideb2 + jj - 1) as usize] = rho * du0[(ideb2 + jj - 1) as usize];
                su1[(ideb3 + jj - 1) as usize] = rho * su1[(ideb3 + jj - 1) as usize];
                du1[(ideb4 + jj - 1) as usize] = rho * du1[(ideb4 + jj - 1) as usize];
                jj += 1;
            }
            iu += 1;
        }

        // tables SomTab and Diftab of discretization of isos V=V0 and V=V1
        // SU0, SU1, DU0 and DU1 are dimensioned in FORTRAN at
        // (1+NBPNTU/2)*NDIMEN*(IORDRV+1)

        size = ((1 + nbpntu / 2) * ndimen) as usize;

        let iso_v0 = constraints.iso_v(self.my_u0, self.my_u1, self.my_v0);
        let mut sv0 = vec![0.0f64; size * (iordrv + 1) as usize];
        sv0.copy_from_slice(&iso_v0.som_tab().expect("Patch::Discretise : SomTab"));
        let iso_v0 = constraints.iso_v(self.my_u0, self.my_u1, self.my_v0);
        let mut dv0 = vec![0.0f64; size * (iordrv + 1) as usize];
        dv0.copy_from_slice(&iso_v0.dif_tab().expect("Patch::Discretise : DifTab"));
        let iso_v1 = constraints.iso_v(self.my_u0, self.my_u1, self.my_v1);
        let mut sv1 = vec![0.0f64; size * (iordrv + 1) as usize];
        sv1.copy_from_slice(&iso_v1.som_tab().expect("Patch::Discretise : SomTab"));
        let iso_v1 = constraints.iso_v(self.my_u0, self.my_u1, self.my_v1);
        let mut dv1 = vec![0.0f64; size * (iordrv + 1) as usize];
        dv1.copy_from_slice(&iso_v1.dif_tab().expect("Patch::Discretise : DifTab"));

        // normalisation
        let mut iv: i32 = 1;
        while iv <= iordrv {
            rho = dv.powi(iv);
            ideb1 = iv * size as i32;
            ideb2 = iv * size as i32;
            ideb3 = iv * size as i32;
            ideb4 = iv * size as i32;
            jj = 1;
            while jj <= size as i32 {
                sv0[(ideb1 + jj - 1) as usize] = rho * sv0[(ideb1 + jj - 1) as usize];
                dv0[(ideb2 + jj - 1) as usize] = rho * dv0[(ideb2 + jj - 1) as usize];
                sv1[(ideb3 + jj - 1) as usize] = rho * sv1[(ideb3 + jj - 1) as usize];
                dv1[(ideb4 + jj - 1) as usize] = rho * dv1[(ideb4 + jj - 1) as usize];
                jj += 1;
            }
            iv += 1;
        }

        // SOSOTB and DIDITB are dimensioned in FORTRAN at
        // (0:NBPNTU/2,0:NBPNTV/2,NDIMEN)

        size = ((1 + nbpntu / 2) * (1 + nbpntv / 2) * ndimen) as usize;

        let mut hsoso = vec![0.0f64; size]; // HSOSO->Init(0.)
        let mut hdidi = vec![0.0f64; size]; // HDIDI->Init(0.)

        // SODITB and DISOTB are dimensioned in FORTRAN at
        // (1:NBPNTU/2,1:NBPNTV/2,NDIMEN)

        size = ((nbpntu / 2) * (nbpntv / 2) * ndimen) as usize;

        let mut hsodi = vec![0.0f64; size]; // HSODI->Init(0.)
        let mut hdiso = vec![0.0f64; size]; // HDISO->Init(0.)

        let mut iercod: i32 = 0;

        //  discretization of polynoms of interpolation
        super::approxf2var_d::mma2cdi_(
            &mut ndimen,
            &mut nbpntu,
            uroot,
            &mut nbpntv,
            vroot,
            &mut iordru,
            &mut iordrv,
            c1,
            c2,
            c3,
            c4,
            &su0,
            &su1,
            &du0,
            &du1,
            &sv0,
            &sv1,
            &dv0,
            &dv1,
            &mut hsoso,
            &mut hsodi,
            &mut hdiso,
            &mut hdidi,
            &mut iercod,
        );

        //  discretization of the square
        let mut udbfn = [0.0f64; 2];
        let mut vdbfn = [0.0f64; 2];
        udbfn[0] = self.my_u0;
        udbfn[1] = self.my_u1;
        vdbfn[0] = self.my_v0;
        vdbfn[1] = self.my_v1;

        size = nbpntu.max(nbpntv) as usize;
        let mut tab = vec![0.0f64; size]; // HTABLE = new ...(1, SIZE)
        let mut pts = vec![0.0f64; size * ndimen as usize]; // HPOINTS

        // GCC 3.0 comment dropped; OCCT passes Func as the evaluator.
        super::approxf2var_f::mma2ds1_(
            &mut ndimen,
            &udbfn,
            &vdbfn,
            func,
            &mut nbpntu,
            &mut nbpntv,
            uroot,
            vroot,
            &mut isofav,
            &mut hsoso,
            &mut hdiso,
            &mut hsodi,
            &mut hdidi,
            &mut pts,
            &mut tab,
            &mut iercod,
        );

        // the results are stored
        if iercod == 0 {
            self.my_disc_is_done = true;
            self.my_soso_tab = Some(hsoso);
            self.my_diso_tab = Some(hdiso);
            self.my_sodi_tab = Some(hsodi);
            self.my_didi_tab = Some(hdidi);
        } else {
            self.my_disc_is_done = false;
        }
    }

    /// OCCT HasResult() (AdvApp2Var_Patch.cxx L374-377).
    pub fn has_result(&self) -> bool {
        self.my_has_result
    }

    /// OCCT IsApproximated() (AdvApp2Var_Patch.cxx L381-384).
    pub fn is_approximated(&self) -> bool {
        self.my_appr_is_done
    }

    /// OCCT AddConstraints(Conditions, Constraints)
    /// (AdvApp2Var_Patch.cxx L388-651).
    pub fn add_constraints(&mut self, conditions: &Context, constraints: &Framework) {
        // data stored in the  Context
        let ndimen: i32 = conditions.total_dimension();
        let mut iercod: i32 = 0;
        let ncflmu = conditions.u_limit();
        let ncflmv = conditions.v_limit();
        let ndeg_u = ncflmu - 1;
        let ndeg_v = ncflmv - 1;

        // data relative to the patch
        let iordru = self.my_ord_in_u;
        let iordrv = self.my_ord_in_v;
        // PATCAN = &myEquation->ChangeArray1()(Lower()) - mutable raw base.
        let patcan = self
            .my_equation
            .as_mut()
            .expect("Patch::AddConstraints : myEquation");

        // curves of approximation of Isos U
        let mut size = (ncflmv * ndimen) as usize;
        // HIsoU0 = ...(1, SIZE * (IORDRU + 1)); contents copied from the iso
        // Polynom (the OCCT ChangeArray1() = Array1() assignment).
        let mut isou0 = vec![0.0f64; size * (iordru + 1) as usize];
        isou0.copy_from_slice(
            &constraints
                .iso_u(self.my_u0, self.my_v0, self.my_v1)
                .polynom()
                .expect("Patch::AddConstraints : Polynom"),
        );
        let mut ncfu0 = vec![0i32; (iordru + 1) as usize];
        // HCFU0->Init(NbCoeff of IsoU(myU0, myV0, myV1)).
        let iso_u0_nbcoeff = constraints.iso_u(self.my_u0, self.my_v0, self.my_v1).nb_coeff();
        ncfu0.fill(iso_u0_nbcoeff);

        let mut isou1 = vec![0.0f64; size * (iordru + 1) as usize];
        isou1.copy_from_slice(
            &constraints
                .iso_u(self.my_u1, self.my_v0, self.my_v1)
                .polynom()
                .expect("Patch::AddConstraints : Polynom"),
        );
        let mut ncfu1 = vec![0i32; (iordru + 1) as usize];
        let iso_u1_nbcoeff = constraints.iso_u(self.my_u1, self.my_v0, self.my_v1).nb_coeff();
        ncfu1.fill(iso_u1_nbcoeff);

        // normalization of Isos U
        let du = (self.my_u1 - self.my_u0) / 2.0;
        let dv = (self.my_v1 - self.my_v0) / 2.0;
        let mut rho: f64;
        let mut valnorm: f64;
        let mut ideb0: i32;
        let mut ideb1: i32;
        let mut jj: i32;

        let mut iu: i32 = 1;
        while iu <= iordru {
            rho = du.powi(iu);
            ideb0 = iu * size as i32;
            ideb1 = iu * size as i32;
            jj = 1;
            while jj <= size as i32 {
                isou0[(ideb0 + jj - 1) as usize] = rho * isou0[(ideb0 + jj - 1) as usize];
                isou1[(ideb1 + jj - 1) as usize] = rho * isou1[(ideb1 + jj - 1) as usize];
                jj += 1;
            }
            iu += 1;
        }

        // curves of approximation of Isos V
        size = (ncflmu * ndimen) as usize;
        let mut isov0 = vec![0.0f64; size * (iordrv + 1) as usize];
        isov0.copy_from_slice(
            &constraints
                .iso_v(self.my_u0, self.my_u1, self.my_v0)
                .polynom()
                .expect("Patch::AddConstraints : Polynom"),
        );
        let mut ncfv0 = vec![0i32; (iordrv + 1) as usize];
        let iso_v0_nbcoeff = constraints.iso_v(self.my_u0, self.my_u1, self.my_v0).nb_coeff();
        ncfv0.fill(iso_v0_nbcoeff);

        let mut isov1 = vec![0.0f64; size * (iordrv + 1) as usize];
        isov1.copy_from_slice(
            &constraints
                .iso_v(self.my_u0, self.my_u1, self.my_v1)
                .polynom()
                .expect("Patch::AddConstraints : Polynom"),
        );
        let mut ncfv1 = vec![0i32; (iordrv + 1) as usize];
        let iso_v1_nbcoeff = constraints.iso_v(self.my_u0, self.my_u1, self.my_v1).nb_coeff();
        ncfv1.fill(iso_v1_nbcoeff);

        //  normalization of Isos V
        let mut iv: i32 = 1;
        while iv <= iordrv {
            rho = dv.powi(iv);
            ideb0 = iv * size as i32;
            ideb1 = iv * size as i32;
            jj = 1;
            while jj <= size as i32 {
                isov0[(ideb0 + jj - 1) as usize] = rho * isov0[(ideb0 + jj - 1) as usize];
                isov1[(ideb1 + jj - 1) as usize] = rho * isov1[(ideb1 + jj - 1) as usize];
                jj += 1;
            }
            iv += 1;
        }

        // add constraints to constant V
        // HHERMV = ...(1, (2*IORDRV+2)*(2*IORDRV+2)).
        let mut hermv = vec![0.0f64; ((2 * iordrv + 2) * (2 * iordrv + 2)) as usize];
        if iordrv >= 0 {
            super::approxf2var_a::mma1her_(&iordrv, &mut hermv, &mut iercod);
            if iercod != 0 {
                panic!("AdvApp2Var_Patch::AddConstraints : Error in FORTRAN");
            }
            super::approxf2var_a::mma2ac2_(
                &ndimen,
                &ndeg_u,
                &ndeg_v,
                &iordrv,
                &ncflmu,
                &ncfv0,
                &isov0,
                &ncfv1,
                &isov1,
                &hermv,
                patcan,
            );
        }

        // add constraints to constant U
        // HHERMU = ...(1, (2*IORDRU+2)*(2*IORDRU+2)).
        let mut hermu = vec![0.0f64; ((2 * iordru + 2) * (2 * iordru + 2)) as usize];
        if iordru >= 0 {
            super::approxf2var_a::mma1her_(&iordru, &mut hermu, &mut iercod);
            if iercod != 0 {
                panic!("AdvApp2Var_Patch::AddConstraints : Error in FORTRAN");
            }
            super::approxf2var_a::mma2ac3_(
                &ndimen,
                &ndeg_u,
                &ndeg_v,
                &iordru,
                &ncflmv,
                &ncfu0,
                &isou0,
                &ncfu1,
                &isou1,
                &hermu,
                patcan,
            );
        }

        // add constraints at the corners
        let mut ideb: i32;
        size = (ndimen * (iordru + 2) * (iordrv + 2)) as usize;
        // HCOINS = ...(1, SIZE * 4): the OCCT array is left uninitialized
        // here (only the (0..=myOrdInU)x(0..=myOrdInV) blocks are filled;
        // never read beyond them for myOrdInU == IORDRU); the rcad buffer
        // zero-fills.
        let mut hcoins = vec![0.0f64; size * 4];

        iu = 0;
        while iu <= self.my_ord_in_u {
            iv = 0;
            while iv <= self.my_ord_in_v {
                rho = du.powi(iu) * dv.powi(iv);

                // -F(U0,V0) and its derivatives normalized on (-1,1)
                ideb = ndimen * iu + ndimen * (iordru + 2) * iv - 1;
                valnorm = -rho * constraints.node(self.my_u0, self.my_v0).point(iu, iv).x;
                hcoins[(1 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u0, self.my_v0).point(iu, iv).y;
                hcoins[(2 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u0, self.my_v0).point(iu, iv).z;
                hcoins[(3 + ideb - 1) as usize] = valnorm;

                // -F(U1,V0) and its derivatives normalized on (-1,1)
                ideb += size as i32;
                valnorm = -rho * constraints.node(self.my_u1, self.my_v0).point(iu, iv).x;
                hcoins[(1 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u1, self.my_v0).point(iu, iv).y;
                hcoins[(2 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u1, self.my_v0).point(iu, iv).z;
                hcoins[(3 + ideb - 1) as usize] = valnorm;

                // -F(U0,V1) and its derivatives normalized on (-1,1)
                ideb += size as i32;
                valnorm = -rho * constraints.node(self.my_u0, self.my_v1).point(iu, iv).x;
                hcoins[(1 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u0, self.my_v1).point(iu, iv).y;
                hcoins[(2 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u0, self.my_v1).point(iu, iv).z;
                hcoins[(3 + ideb - 1) as usize] = valnorm;

                // -F(U1,V1) and its derivatives normalized on (-1,1)
                ideb += size as i32;
                valnorm = -rho * constraints.node(self.my_u1, self.my_v1).point(iu, iv).x;
                hcoins[(1 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u1, self.my_v1).point(iu, iv).y;
                hcoins[(2 + ideb - 1) as usize] = valnorm;
                valnorm = -rho * constraints.node(self.my_u1, self.my_v1).point(iu, iv).z;
                hcoins[(3 + ideb - 1) as usize] = valnorm;
                iv += 1;
            }
            iu += 1;
        }

        //  tables required for FORTRAN
        let iordmx = iordru.max(iordrv);
        let mut extr = vec![0.0f64; (2 * iordmx + 2) as usize];
        let mut fact = vec![0.0f64; (iordmx + 1) as usize];

        let mut idim: i32;
        let mut ncf0: i32;
        let mut ncf1: i32;
        let mut iun: i32 = 1;
        let mut is_off: usize;

        // add extremities of isos U
        iu = 1;
        while iu <= iordru + 1 {
            ncf0 = ncfu0[(iu - 1) as usize];
            ncf1 = ncfu1[(iu - 1) as usize];
            idim = 1;
            while idim <= ndimen {
                // OCCT: Is = IsoU0 + NCFLMV*(idim-1) + NCFLMV*NDIMEN*(iu-1).
                is_off = (ncflmv * (idim - 1) + ncflmv * ndimen * (iu - 1)) as usize;
                super::math_base::mmdrc11_(
                    &mut { iordrv },
                    &mut iun,
                    &mut ncf0,
                    &isou0[is_off..],
                    &mut extr,
                    &mut fact,
                );
                iv = 1;
                while iv <= iordrv + 1 {
                    // OCCT: ideb = Lower + NDIMEN*(iu-1) +
                    //   NDIMEN*(IORDRU+2)*(iv-1) - 1.
                    ideb = ndimen * (iu - 1) + ndimen * (iordru + 2) * (iv - 1) - 1;
                    hcoins[(idim + ideb) as usize] +=
                        extr[(1 + 2 * (iv - 1) - 1) as usize];
                    hcoins[(2 * size as i32 + idim + ideb) as usize] +=
                        extr[(2 + 2 * (iv - 1) - 1) as usize];
                    iv += 1;
                }
                is_off = (ncflmv * (idim - 1) + ncflmv * ndimen * (iu - 1)) as usize;
                super::math_base::mmdrc11_(
                    &mut { iordrv },
                    &mut iun,
                    &mut ncf1,
                    &isou1[is_off..],
                    &mut extr,
                    &mut fact,
                );
                iv = 1;
                while iv <= iordrv + 1 {
                    ideb = ndimen * (iu - 1) + ndimen * (iordru + 2) * (iv - 1) - 1;
                    hcoins[(size as i32 + idim + ideb) as usize] +=
                        extr[(1 + 2 * (iv - 1) - 1) as usize];
                    hcoins[(3 * size as i32 + idim + ideb) as usize] +=
                        extr[(2 + 2 * (iv - 1) - 1) as usize];
                    iv += 1;
                }
                idim += 1;
            }
            iu += 1;
        }

        // add extremities of isos V
        iv = 1;
        while iv <= iordrv + 1 {
            ncf0 = ncfv0[(iv - 1) as usize];
            ncf1 = ncfv1[(iv - 1) as usize];
            idim = 1;
            while idim <= ndimen {
                // OCCT: Is = IsoV0 + NCFLMU*(idim-1) + NCFLMU*NDIMEN*(iv-1).
                is_off = (ncflmu * (idim - 1) + ncflmu * ndimen * (iv - 1)) as usize;
                super::math_base::mmdrc11_(
                    &mut { iordru },
                    &mut iun,
                    &mut ncf0,
                    &isov0[is_off..],
                    &mut extr,
                    &mut fact,
                );
                iu = 1;
                while iu <= iordru + 1 {
                    ideb = ndimen * (iu - 1) + ndimen * (iordru + 2) * (iv - 1) - 1;
                    hcoins[(idim + ideb) as usize] += extr[(1 + 2 * (iu - 1) - 1) as usize];
                    hcoins[(size as i32 + idim + ideb) as usize] +=
                        extr[(2 + 2 * (iu - 1) - 1) as usize];
                    iu += 1;
                }
                is_off = (ncflmu * (idim - 1) + ncflmu * ndimen * (iv - 1)) as usize;
                super::math_base::mmdrc11_(
                    &mut { iordru },
                    &mut iun,
                    &mut ncf1,
                    &isov1[is_off..],
                    &mut extr,
                    &mut fact,
                );
                iu = 1;
                while iu <= iordru + 1 {
                    ideb = ndimen * (iu - 1) + ndimen * (iordru + 2) * (iv - 1) - 1;
                    hcoins[(2 * size as i32 + idim + ideb) as usize] +=
                        extr[(1 + 2 * (iu - 1) - 1) as usize];
                    hcoins[(3 * size as i32 + idim + ideb) as usize] +=
                        extr[(2 + 2 * (iu - 1) - 1) as usize];
                    iu += 1;
                }
                idim += 1;
            }
            iv += 1;
        }

        // add all to PATCAN
        // C1 = base; C2 = C1 + SIZE; C3 = C2 + SIZE; C4 = C3 + SIZE.
        let c1: &[f64] = &hcoins[0..size];
        let c2: &[f64] = &hcoins[size..2 * size];
        let c3: &[f64] = &hcoins[2 * size..3 * size];
        let c4: &[f64] = &hcoins[3 * size..4 * size];
        if iordru >= 0 && iordrv >= 0 {
            super::approxf2var_a::mma2ac1_(
                &ndimen,
                &ndeg_u,
                &ndeg_v,
                &iordru,
                &iordrv,
                c1,
                c2,
                c3,
                c4,
                &hermu,
                &hermv,
                patcan,
            );
        }
    }

    /// OCCT AddErrors(Constraints) (AdvApp2Var_Patch.cxx L655-760).
    pub fn add_errors(&mut self, constraints: &Framework) {
        let nbsesp: i32 = 1;
        let mut iesp: i32;
        let mut iu: i32;
        let mut iv: i32;

        let mut err_u: f64;
        let mut err_v: f64;
        let mut error: f64;
        let hmax: [f64; 4] = [0.0, 1.0, 1.5, 1.75];

        // myMaxErrors / myMoyErrors / myIsoErrors handle dereferences (set
        // by MakeApprox before AddErrors is called).
        iesp = 1;
        while iesp <= nbsesp {
            //  error max in sub-space iesp
            err_u = 0.;
            iv = 1;
            while iv <= self.my_ord_in_v + 1 {
                error = constraints
                    .iso_v(self.my_u0, self.my_u1, self.my_v0)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, iv);
                err_u = err_u.max(error);
                error = constraints
                    .iso_v(self.my_u0, self.my_u1, self.my_v1)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, iv);
                err_u = err_u.max(error);
                iv += 1;
            }
            err_v = 0.;
            iu = 1;
            while iu <= self.my_ord_in_u + 1 {
                error = constraints
                    .iso_u(self.my_u0, self.my_v0, self.my_v1)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, iu);
                err_v = err_v.max(error);
                error = constraints
                    .iso_u(self.my_u1, self.my_v0, self.my_v1)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, iu);
                err_v = err_v.max(error);
                iu += 1;
            }
            let my_max_errors = self.my_max_errors.as_mut().expect("Patch::AddErrors");
            my_max_errors[(iesp - 1) as usize] +=
                err_u * hmax[(self.my_ord_in_v + 1) as usize] + err_v * hmax[(self.my_ord_in_u + 1) as usize];

            // average error in sub-space iesp
            err_u = 0.;
            iv = 1;
            while iv <= self.my_ord_in_v + 1 {
                error = constraints
                    .iso_v(self.my_u0, self.my_u1, self.my_v0)
                    .moy_errors()
                    .expect("Patch::AddErrors : MoyErrors")
                    .value(iesp, iv);
                err_u = err_u.max(error);
                error = constraints
                    .iso_v(self.my_u0, self.my_u1, self.my_v1)
                    .moy_errors()
                    .expect("Patch::AddErrors : MoyErrors")
                    .value(iesp, iv);
                err_u = err_u.max(error);
                iv += 1;
            }
            err_v = 0.;
            iu = 1;
            while iu <= self.my_ord_in_u + 1 {
                error = constraints
                    .iso_u(self.my_u0, self.my_v0, self.my_v1)
                    .moy_errors()
                    .expect("Patch::AddErrors : MoyErrors")
                    .value(iesp, iu);
                err_v = err_v.max(error);
                error = constraints
                    .iso_u(self.my_u1, self.my_v0, self.my_v1)
                    .moy_errors()
                    .expect("Patch::AddErrors : MoyErrors")
                    .value(iesp, iu);
                err_v = err_v.max(error);
                iu += 1;
            }
            let my_moy_errors = self.my_moy_errors.as_mut().expect("Patch::AddErrors");
            error = my_moy_errors[(iesp - 1) as usize];
            error *= error;
            error += err_u * hmax[(self.my_ord_in_v + 1) as usize] * err_u * hmax[(self.my_ord_in_v + 1) as usize]
                + err_v * hmax[(self.my_ord_in_u + 1) as usize] * err_v * hmax[(self.my_ord_in_u + 1) as usize];
            my_moy_errors[(iesp - 1) as usize] = error.sqrt();

            // max errors at iso-borders
            // HERISO = new NCollection_HArray2<double>(1, NBSESP, 1, 4).
            let mut heriso = Array2::new(1, nbsesp, 1, 4);
            heriso.set_value(
                iesp,
                1,
                constraints
                    .iso_v(self.my_u0, self.my_u1, self.my_v0)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, 1),
            );
            heriso.set_value(
                iesp,
                2,
                constraints
                    .iso_v(self.my_u0, self.my_u1, self.my_v1)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, 1),
            );
            heriso.set_value(
                iesp,
                3,
                constraints
                    .iso_u(self.my_u0, self.my_v0, self.my_v1)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, 1),
            );
            heriso.set_value(
                iesp,
                4,
                constraints
                    .iso_u(self.my_u1, self.my_v0, self.my_v1)
                    .max_errors()
                    .expect("Patch::AddErrors : MaxErrors")
                    .value(iesp, 1),
            );

            // calculate max errors at the corners
            let mut emax1: f64 = 0.;
            let mut emax2: f64 = 0.;
            let mut emax3: f64 = 0.;
            let mut emax4: f64 = 0.;
            iu = 0;
            while iu <= self.my_ord_in_u {
                iv = 0;
                while iv <= self.my_ord_in_v {
                    error = constraints.node(self.my_u0, self.my_v0).error(iu, iv);
                    emax1 = emax1.max(error);
                    error = constraints.node(self.my_u1, self.my_v0).error(iu, iv);
                    emax2 = emax2.max(error);
                    error = constraints.node(self.my_u0, self.my_v1).error(iu, iv);
                    emax3 = emax3.max(error);
                    error = constraints.node(self.my_u1, self.my_v1).error(iu, iv);
                    emax4 = emax4.max(error);
                    iv += 1;
                }
                iu += 1;
            }

            // calculate max errors on borders
            let err1 = emax1.max(emax2);
            let err2 = emax3.max(emax4);
            let err3 = emax1.max(emax3);
            let err4 = emax2.max(emax4);

            //   calculate final errors on internal isos
            if constraints.iso_v(self.my_u0, self.my_u1, self.my_v0).position() == 0 {
                *heriso.change_value(iesp, 1) += err1 * hmax[(self.my_ord_in_u + 1) as usize];
            }
            if constraints.iso_v(self.my_u0, self.my_u1, self.my_v1).position() == 0 {
                *heriso.change_value(iesp, 2) += err2 * hmax[(self.my_ord_in_u + 1) as usize];
            }
            if constraints.iso_u(self.my_u0, self.my_v0, self.my_v1).position() == 0 {
                *heriso.change_value(iesp, 3) += err3 * hmax[(self.my_ord_in_v + 1) as usize];
            }
            if constraints.iso_u(self.my_u1, self.my_v0, self.my_v1).position() == 0 {
                *heriso.change_value(iesp, 4) += err4 * hmax[(self.my_ord_in_v + 1) as usize];
            }
            self.my_iso_errors = Some(heriso);
            iesp += 1;
        }
    }

    /// OCCT MakeApprox(Conditions, Constraints, NumDec)
    /// (AdvApp2Var_Patch.cxx L764-964).  The dead initializers of
    /// NDMINU/NDMINV (= 1, cxx L796) and ITYDEC/IERCOD/NCOEFU/NCOEFV
    /// (= 0, cxx L851-854) are OCCT-faithful.
    #[allow(unused_assignments)]
    pub fn make_approx(
        &mut self,
        conditions: &Context,
        constraints: &Framework,
        num_dec: i32,
    ) {
        // data stored in the Context
        let mut ndimen: i32;
        let nbsesp: i32;
        let ndimse: i32;
        let mut nbpntu: i32;
        let mut nbpntv: i32;
        let ncflmu: i32;
        let ncflmv: i32;
        let ndjacu: i32;
        let ndjacv: i32;
        let ndeg_u: i32;
        let ndeg_v: i32;
        let njacu: i32;
        let njacv: i32;
        ndimen = conditions.total_dimension();
        nbsesp = conditions.total_number_ssp();
        ndimse = 3;
        nbpntu = conditions.u_roots().expect("Patch::MakeApprox : URoots").len() as i32;
        if self.my_ord_in_u > -1 {
            nbpntu -= 2;
        }
        nbpntv = conditions.v_roots().expect("Patch::MakeApprox : VRoots").len() as i32;
        if self.my_ord_in_v > -1 {
            nbpntv -= 2;
        }
        ncflmu = conditions.u_limit();
        ncflmv = conditions.v_limit();
        ndeg_u = ncflmu - 1;
        ndeg_v = ncflmv - 1;
        ndjacu = conditions.u_jac_deg();
        ndjacv = conditions.v_jac_deg();
        njacu = ndjacu + 1;
        njacv = ndjacv + 1;

        // data relative to the processed patch
        let mut iordru = self.my_ord_in_u;
        let mut iordrv = self.my_ord_in_v;
        // NDMINU and NDMINV depend on the nb of coeff of neighboring isos
        // and of the required order of continuity
        let mut ndminu: i32 = 1;
        let mut ndminv: i32 = 1;
        let mut ncoefu: i32;
        let mut ncoefv: i32;
        ndminu = 1.max(2 * iordru + 1);
        ncoefu = constraints.iso_v(self.my_u0, self.my_u1, self.my_v0).nb_coeff() - 1;
        ndminu = ndminu.max(ncoefu);
        ncoefu = constraints.iso_v(self.my_u0, self.my_u1, self.my_v1).nb_coeff() - 1;
        ndminu = ndminu.max(ncoefu);

        ndminv = 1.max(2 * iordrv + 1);
        ncoefv = constraints.iso_u(self.my_u0, self.my_v0, self.my_v1).nb_coeff() - 1;
        ndminv = ndminv.max(ncoefv);
        ncoefv = constraints.iso_u(self.my_u1, self.my_v0, self.my_v1).nb_coeff() - 1;
        ndminv = ndminv.max(ncoefv);

        // tables of approximations
        // HEPSAPR = ...(1, NBSESP); HEPSFRO = ...(1, NBSESP * 8).
        let mut epsapr = vec![0.0f64; nbsesp as usize];
        let mut epsfro = vec![0.0f64; (nbsesp * 8) as usize];
        let itoler = conditions.i_toler().expect("Patch::MakeApprox : IToler");
        let ftoler = conditions.f_toler().expect("Patch::MakeApprox : FToler");
        let ctoler = conditions.c_toler().expect("Patch::MakeApprox : CToler");
        let mut iesp: i32 = 1;
        while iesp <= nbsesp {
            epsapr[(iesp - 1) as usize] = itoler[(iesp - 1) as usize];
            epsfro[(iesp - 1) as usize] = ftoler.value(iesp, 1);
            epsfro[(iesp + nbsesp - 1) as usize] = ftoler.value(iesp, 2);
            epsfro[(iesp + 2 * nbsesp - 1) as usize] = ftoler.value(iesp, 3);
            epsfro[(iesp + 3 * nbsesp - 1) as usize] = ftoler.value(iesp, 4);
            epsfro[(iesp + 4 * nbsesp - 1) as usize] = ctoler.value(iesp, 1);
            epsfro[(iesp + 5 * nbsesp - 1) as usize] = ctoler.value(iesp, 2);
            epsfro[(iesp + 6 * nbsesp - 1) as usize] = ctoler.value(iesp, 3);
            epsfro[(iesp + 7 * nbsesp - 1) as usize] = ctoler.value(iesp, 4);
            iesp += 1;
        }

        let mut size = ((1 + ndjacu) * (1 + ndjacv) * ndimen) as usize;
        // HPJAC = ...(1, SIZE); PATAUX = ...(1, 2 * SIZE).
        let mut patjac = vec![0.0f64; size];
        let mut pataux = vec![0.0f64; 2 * size];
        size = (ncflmu * ncflmv * ndimen) as usize;
        let mut patcan = vec![0.0f64; size];
        let mut herrmax = vec![0.0f64; nbsesp as usize];
        let mut herrmoy = vec![0.0f64; nbsesp as usize];

        // tables of discretization of the square (the Discretise results;
        // the OCCT bodies write through the stored arrays - the rcad copies
        // are unobservable-equal since only Discretise/MakeApprox consume
        // these tables).
        let mut sosotb = self.my_soso_tab.clone().expect("Patch::MakeApprox : SosoTab");
        let mut disotb = self.my_diso_tab.clone().expect("Patch::MakeApprox : DisoTab");
        let mut soditb = self.my_sodi_tab.clone().expect("Patch::MakeApprox : SodiTab");
        let mut diditb = self.my_didi_tab.clone().expect("Patch::MakeApprox : DidiTab");

        //  approximation
        let mut itydec: i32 = 0;
        let mut iercod: i32 = 0;
        let iun: i32 = 1;
        let itrois: i32 = 3;
        ncoefu = 0;
        ncoefv = 0;
        // OCCT: mma2ce1_((int*)&NumDec, ...) - the const parameter is passed
        // by address; rcad copies it into a mutable local.
        let mut numdec = num_dec;
        super::approxf2var_e::mma2ce1_(
            &mut numdec,
            &mut ndimen,
            &mut { nbsesp },
            std::slice::from_ref(&ndimse),
            &mut ndminu,
            &mut ndminv,
            &mut { ndeg_u },
            &mut { ndeg_v },
            &mut { ndjacu },
            &mut { ndjacv },
            &mut iordru,
            &mut iordrv,
            &mut nbpntu,
            &mut nbpntv,
            &epsapr,
            &mut sosotb,
            &mut disotb,
            &mut soditb,
            &mut diditb,
            &mut patjac,
            &mut herrmax,
            &mut herrmoy,
            &mut ncoefu,
            &mut ncoefv,
            &mut itydec,
            &mut iercod,
        );

        // results
        self.my_cut_sense = itydec;
        if itydec == 0 && iercod <= 0 {
            self.my_has_result = true;
            self.my_appr_is_done = iercod == 0;
            self.my_nb_coeff_in_u = ncoefu + 1;
            self.my_nb_coeff_in_v = ncoefv + 1;
            self.my_max_errors = Some(herrmax);
            self.my_moy_errors = Some(herrmoy);

            // Passage to canonic on [-1,1]
            // OCCT passes PATJAC twice (aliased in/out buffer); the mmfmca9_
            // write index never exceeds its read index (the target dims
            // NCOEFU/NCOEFV are <= the source dims NJacU/NJacV), so the
            // in-place transform equals the snapshot semantics reproduced
            // here for the & / &mut split.
            let patjac_snap = patjac.clone();
            super::math_base_b::mmfmca9_(
                &njacu,
                &njacv,
                &ndimen,
                &self.my_nb_coeff_in_u,
                &self.my_nb_coeff_in_v,
                &ndimen,
                &patjac_snap,
                &mut patjac,
            );
            super::approxf2var_d::mma2can_(
                &ncflmu,
                &ncflmv,
                &ndimen,
                &self.my_ord_in_u,
                &self.my_ord_in_v,
                &self.my_nb_coeff_in_u,
                &self.my_nb_coeff_in_v,
                &patjac,
                &mut pataux,
                &mut patcan,
                &mut iercod,
            );
            if iercod != 0 {
                panic!("AdvApp2Var_Patch::MakeApprox : Error in FORTRAN");
            }
            self.my_equation = Some(patcan);

            // Add constraints and errors
            self.add_constraints(conditions, constraints);
            self.add_errors(constraints);

            // Reduction of degrees if possible
            // PATCAN = &myEquation->ChangeArray1()(Lower()) (re-fetched).
            // OCCT passes &myNbCoeffInU / &myNbCoeffInV by address; the
            // values are copied in and out around the rcad call.  The two
            // &iun arguments (nbupat / nbvpat) alias the same OCCT
            // variable; rcad passes two fresh copies of it.  ERRMAX is the
            // buffer stored in myMaxErrors (the OCCT handle is shared), so
            // the degree-reduction update lands in the member directly.
            let patcan = self
                .my_equation
                .as_mut()
                .expect("Patch::MakeApprox : myEquation");
            let errmax_ref = self.my_max_errors.as_mut().expect("Patch::MakeApprox");
            let mut nbcoeffu = self.my_nb_coeff_in_u;
            let mut nbcoeffv = self.my_nb_coeff_in_v;
            super::approxf2var_b::mma2fx6_(
                &mut { ncflmu },
                &mut { ncflmv },
                &mut { ndimen },
                &mut { nbsesp },
                std::slice::from_ref(&itrois),
                &mut { iun },
                &mut { iun },
                &mut iordru,
                &mut iordrv,
                &epsapr,
                &epsfro,
                patcan,
                errmax_ref,
                std::slice::from_mut(&mut nbcoeffu),
                std::slice::from_mut(&mut nbcoeffv),
            );
            self.my_nb_coeff_in_u = nbcoeffu;
            self.my_nb_coeff_in_v = nbcoeffv;

            // transposition (NCFLMU,NCFLMV,NDIMEN)Fortran-C++
            let mut a_iu: i32;
            let mut a_in: i32;
            let equation = self.my_equation.as_ref().expect("Patch::MakeApprox");
            let mut hpaux_buf = pataux;
            let mut dim: i32 = 1;
            while dim <= ndimen {
                a_in = (dim - 1) * ncflmu * ncflmv;
                let mut ii: i32 = 1;
                while ii <= ncflmu {
                    a_iu = (ii - 1) * ndimen * ncflmv;
                    let mut jj: i32 = 1;
                    while jj <= ncflmv {
                        // OCCT: HPAUX->SetValue(dim + NDIMEN*(jj-1) + aIU,
                        //   myEquation->Value(ii + NCFLMU*(jj-1) + aIN)).
                        hpaux_buf[(dim + ndimen * (jj - 1) + a_iu - 1) as usize] =
                            equation[(ii + ncflmu * (jj - 1) + a_in - 1) as usize];
                        jj += 1;
                    }
                    ii += 1;
                }
                dim += 1;
            }
            self.my_equation = Some(hpaux_buf);
        } else {
            self.my_appr_is_done = false;
            self.my_has_result = false;
        }
    }

    /// OCCT ChangeDomain(a, b, c, d) (AdvApp2Var_Patch.cxx L968-974).
    pub fn change_domain(&mut self, a: f64, b: f64, c: f64, d: f64) {
        self.my_u0 = a;
        self.my_u1 = b;
        self.my_v0 = c;
        self.my_v1 = d;
    }

    /// OCCT ResetApprox() (AdvApp2Var_Patch.cxx L978-982).
    pub fn reset_approx(&mut self) {
        self.my_appr_is_done = false;
        self.my_has_result = false;
    }

    /// OCCT OverwriteApprox() (AdvApp2Var_Patch.cxx L986-992).
    pub fn overwrite_approx(&mut self) {
        if self.my_has_result {
            self.my_appr_is_done = true;
        }
    }

    /// OCCT U0() (AdvApp2Var_Patch.cxx L996-999).
    pub fn u0(&self) -> f64 {
        self.my_u0
    }

    /// OCCT U1() (AdvApp2Var_Patch.cxx L1003-1006).
    pub fn u1(&self) -> f64 {
        self.my_u1
    }

    /// OCCT V0() (AdvApp2Var_Patch.cxx L1010-1013).
    pub fn v0(&self) -> f64 {
        self.my_v0
    }

    /// OCCT V1() (AdvApp2Var_Patch.cxx L1017-1020).
    pub fn v1(&self) -> f64 {
        self.my_v1
    }

    /// OCCT UOrder() (AdvApp2Var_Patch.cxx L1024-1027).
    pub fn u_order(&self) -> i32 {
        self.my_ord_in_u
    }

    /// OCCT VOrder() (AdvApp2Var_Patch.cxx L1031-1034).
    pub fn v_order(&self) -> i32 {
        self.my_ord_in_v
    }

    /// OCCT CutSense() (AdvApp2Var_Patch.cxx L1038-1041).
    pub fn cut_sense(&self) -> i32 {
        self.my_cut_sense
    }

    /// OCCT CutSense(Crit, NumDec) (AdvApp2Var_Patch.cxx L1045-1063) - the
    /// overload receives a distinct rcad name (Rust has no overloading).
    pub fn cut_sense_criterion(&self, crit: &dyn Criterion, num_dec: i32) -> i32 {
        let crit_rel = crit.crit_type() == CriterionType::Relative;
        if crit_rel && !self.is_approximated() {
            self.my_cut_sense
        } else {
            if crit.is_satisfied(self) {
                0
            } else {
                num_dec
            }
        }
    }

    /// OCCT NbCoeffInU() (AdvApp2Var_Patch.cxx L1067-1070).
    pub fn nb_coeff_in_u(&self) -> i32 {
        self.my_nb_coeff_in_u
    }

    /// OCCT NbCoeffInV() (AdvApp2Var_Patch.cxx L1074-1077).
    pub fn nb_coeff_in_v(&self) -> i32 {
        self.my_nb_coeff_in_v
    }

    /// OCCT ChangeNbCoeff(NbCoeffU, NbCoeffV) (AdvApp2Var_Patch.cxx
    /// L1081-1091).
    pub fn change_nb_coeff(&mut self, nb_coeff_u: i32, nb_coeff_v: i32) {
        if self.my_nb_coeff_in_u < nb_coeff_u {
            self.my_nb_coeff_in_u = nb_coeff_u;
        }
        if self.my_nb_coeff_in_v < nb_coeff_v {
            self.my_nb_coeff_in_v = nb_coeff_v;
        }
    }

    /// OCCT MaxErrors() (AdvApp2Var_Patch.cxx L1095-1098) - handle copy
    /// becomes a clone (None = null handle).
    pub fn max_errors(&self) -> Option<Vec<f64>> {
        self.my_max_errors.clone()
    }

    /// OCCT AverageErrors() (AdvApp2Var_Patch.cxx L1102-1105).
    pub fn average_errors(&self) -> Option<Vec<f64>> {
        self.my_moy_errors.clone()
    }

    /// OCCT IsoErrors() (AdvApp2Var_Patch.cxx L1109-1112).
    pub fn iso_errors(&self) -> Option<Array2<f64>> {
        self.my_iso_errors.clone()
    }

    /// OCCT Poles(SSPIndex, Cond) (AdvApp2Var_Patch.cxx L1116-1145).
    pub fn poles(&self, ssp_index: i32, cond: &Context) -> Array2<DVec3> {
        // SousEquation = myEquation when the subspace selection is trivial.
        if !(cond.total_number_ssp() == 1 && ssp_index == 1) {
            panic!("AdvApp2Var_Patch::Poles :  SSPIndex out of range");
        }
        let sous_equation = self.my_equation.as_ref().expect("Patch::Poles : myEquation");

        // NCollection_Array1<double> Intervalle(1, 2) = (-1, 1).
        let mut intervalle = vec![0.0f64; 2];
        intervalle[0] = -1.0;
        intervalle[1] = 1.0;

        // NCollection_Array1<int> NbCoeff(1, 2).
        let mut nb_coeff = vec![0i32; 2];
        nb_coeff[0] = self.my_nb_coeff_in_u;
        nb_coeff[1] = self.my_nb_coeff_in_v;

        // GAP: Convert_GridPolynomialToPoles (TKMath/Convert,
        // Convert_GridPolynomialToPoles.cxx, the 12-argument grid
        // constructor consumed here and in
        // AdvApp2Var_ApproxAFunc2Var::ConvertBS) is not yet translated - it
        // is owned by the convert_comp_polynomial_to_poles.rs batch.  OCCT
        // anchor: AdvApp2Var_Patch.cxx L1137-1142
        // (Convert_GridPolynomialToPoles Conv(Cond.ULimit()-1,
        // Cond.VLimit()-1, NbCoeff, SousEquation->Array1(), Intervalle,
        // Intervalle); return new HArray2<gp_Pnt>(Conv.Poles())).
        let _ = (&nb_coeff, &sous_equation, &intervalle);
        panic!(
            "GAP: Convert_GridPolynomialToPoles (12-arg grid ctor, \
             AdvApp2Var_Patch.cxx L1137-1144) is not translated yet"
        );
    }

    /// OCCT Coefficients(SSPIndex, Cond) (AdvApp2Var_Patch.cxx L1149-1163).
    pub fn coefficients(&self, ssp_index: i32, cond: &Context) -> Vec<f64> {
        if !(cond.total_number_ssp() == 1 && ssp_index == 1) {
            panic!("AdvApp2Var_Patch::Poles :  SSPIndex out of range");
        }
        // return SousEquation (handle copy -> clone of the equation).
        self.my_equation.clone().expect("Patch::Coefficients : myEquation")
    }

    /// OCCT CritValue() (AdvApp2Var_Patch.cxx L1167-1170).
    pub fn crit_value(&self) -> f64 {
        self.my_crit_value
    }

    /// OCCT SetCritValue(dist) (AdvApp2Var_Patch.cxx L1174-1177).
    pub fn set_crit_value(&mut self, dist: f64) {
        self.my_crit_value = dist;
    }
}
