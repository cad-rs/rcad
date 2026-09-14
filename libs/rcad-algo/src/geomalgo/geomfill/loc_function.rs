//! OCCT GeomFill_LocFunction (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_LocFunction.hxx (L30-63) + GeomFill_LocFunction.cxx (whole file
//! L25-147).
//!
//! Consumed by `GeomFill_Sweep::BuildProduct` (cxx L392-477) through the
//! `GeomFill_Sweep_Eval` evaluator ([`super::SweepEval`]): the function
//! packs the location law frame into four gp_Vec slots — V(1) is the law
//! translation, V(2..4) are the columns of the transform matrix (and
//! likewise DV/D2V for the derivatives) — so the 4-section 3d
//! approximation of BuildProduct carries the formal tensor product.
//!
//! Architecture differences:
//! - `handle(GeomFill_LocationLaw) myLaw` (hxx L56) maps to
//!   `Rc<RefCell<dyn LocationLaw>>` (the law handle shared with the Sweep
//!   consumer).
//! - `NCollection_Array1<gp_Vec>` with bounds (1..4) maps to `[DVec3; 4]`
//!   (0-based: OCCT index i is `self.[i - 1]`).
//! - the .cxx methods work on local `gp_Mat aM / aDM / aD2M`; the hxx
//!   members `M / DM / D2M` (hxx L60-62) are unused there as well and stay
//!   for the 1:1 header mapping.

use std::cell::RefCell;
use std::rc::Rc;

use glam::{DVec2, DVec3};

use super::super::gp_mat::GpMat;
use super::super::location_law::LocationLaw;

/// OCCT GeomFill_LocFunction (hxx L30-63).
pub struct LocFunction {
    /// OCCT handle(GeomFill_LocationLaw) myLaw (hxx L56).
    my_law: Rc<RefCell<dyn LocationLaw>>,
    /// OCCT NCollection_Array1<gp_Vec> V (hxx L57, bounds 1..4).
    v: [DVec3; 4],
    /// OCCT NCollection_Array1<gp_Vec> DV (hxx L58, bounds 1..4).
    dv: [DVec3; 4],
    /// OCCT NCollection_Array1<gp_Vec> D2V (hxx L59, bounds 1..4).
    d2v: [DVec3; 4],
    /// OCCT gp_Mat M (hxx L60) — not referenced by the .cxx.
    #[allow(dead_code)]
    m: GpMat,
    /// OCCT gp_Mat DM (hxx L61) — not referenced by the .cxx.
    #[allow(dead_code)]
    dm: GpMat,
    /// OCCT gp_Mat D2M (hxx L62) — not referenced by the .cxx.
    #[allow(dead_code)]
    d2m: GpMat,
}

impl LocFunction {
    /// OCCT GeomFill_LocFunction::GeomFill_LocFunction (cxx L25-32) —
    /// V(1, 4), DV(1, 4), D2V(1, 4); myLaw = Law.
    pub fn new(law: Rc<RefCell<dyn LocationLaw>>) -> Self {
        LocFunction {
            my_law: law,
            v: [DVec3::ZERO; 4],
            dv: [DVec3::ZERO; 4],
            d2v: [DVec3::ZERO; 4],
            m: GpMat::identity(),
            dm: GpMat::identity(),
            d2m: GpMat::identity(),
        }
    }

    /// OCCT GeomFill_LocFunction::D0 (cxx L34-47) — compute the section for
    /// v = param (the First/Last parameters are commented out upstream).
    pub fn d0(&mut self, param: f64, _first: f64, _last: f64) -> bool {
        let mut a_m = GpMat::identity();
        let b = self.my_law.borrow().d0(param, &mut a_m, &mut self.v[0]);
        self.v[1] = super::gp_mat_column(&a_m, 1);
        self.v[2] = super::gp_mat_column(&a_m, 2);
        self.v[3] = super::gp_mat_column(&a_m, 3);
        b
    }

    /// OCCT GeomFill_LocFunction::D1 (cxx L49-69) — compute the first
    /// derivative in v direction of the section for v = param.
    #[allow(clippy::too_many_arguments)]
    pub fn d1(&mut self, param: f64, _first: f64, _last: f64) -> bool {
        // OCCT: NCollection_Array1<gp_Pnt2d> T1(1, 1);
        //       NCollection_Array1<gp_Vec2d> T2(1, 1).
        let mut t1 = [DVec2::ZERO; 1];
        let mut t2 = [DVec2::ZERO; 1];
        let mut a_m = GpMat::identity();
        let mut a_dm = GpMat::identity();
        let b = self.my_law.borrow().d1(
            param,
            &mut a_m,
            &mut self.v[0],
            &mut a_dm,
            &mut self.dv[0],
            &mut t1,
            &mut t2,
        );

        self.v[1] = super::gp_mat_column(&a_m, 1);
        self.v[2] = super::gp_mat_column(&a_m, 2);
        self.v[3] = super::gp_mat_column(&a_m, 3);

        self.dv[1] = super::gp_mat_column(&a_dm, 1);
        self.dv[2] = super::gp_mat_column(&a_dm, 2);
        self.dv[3] = super::gp_mat_column(&a_dm, 3);
        b
    }

    /// OCCT GeomFill_LocFunction::D2 (cxx L71-104) — compute the second
    /// derivative in v direction of the section for v = param.
    #[allow(clippy::too_many_arguments)]
    pub fn d2(&mut self, param: f64, _first: f64, _last: f64) -> bool {
        // OCCT: NCollection_Array1<gp_Pnt2d> T1(1, 1);
        //       NCollection_Array1<gp_Vec2d> T2(1, 1), T3(1, 1).
        let mut t1 = [DVec2::ZERO; 1];
        let mut t2 = [DVec2::ZERO; 1];
        let mut t3 = [DVec2::ZERO; 1];
        let mut a_m = GpMat::identity();
        let mut a_dm = GpMat::identity();
        let mut a_d2m = GpMat::identity();
        let b = self.my_law.borrow().d2(
            param,
            &mut a_m,
            &mut self.v[0],
            &mut a_dm,
            &mut self.dv[0],
            &mut a_d2m,
            &mut self.d2v[0],
            &mut t1,
            &mut t2,
            &mut t3,
        );

        self.v[1] = super::gp_mat_column(&a_m, 1);
        self.v[2] = super::gp_mat_column(&a_m, 2);
        self.v[3] = super::gp_mat_column(&a_m, 3);

        self.dv[1] = super::gp_mat_column(&a_dm, 1);
        self.dv[2] = super::gp_mat_column(&a_dm, 2);
        self.dv[3] = super::gp_mat_column(&a_dm, 3);

        self.d2v[1] = super::gp_mat_column(&a_d2m, 1);
        self.d2v[2] = super::gp_mat_column(&a_d2m, 2);
        self.d2v[3] = super::gp_mat_column(&a_d2m, 3);

        b
    }

    /// OCCT GeomFill_LocFunction::DN (cxx L106-147) — dispatch to D0/D1/D2
    /// by Order and copy the packed 12 doubles (4 gp_Vec slots) into
    /// Result; Ier = Order + 1 flags a failed (or unsupported) order.
    pub fn dn(
        &mut self,
        param: f64,
        first: f64,
        last: f64,
        order: i32,
        result: &mut [f64],
        ier: &mut i32,
    ) {
        // OCCT: bool B; double* AddrResult = &Result;
        //       const double* LocalResult = nullptr;
        let b;
        let local_result: Option<&[DVec3; 4]>;
        *ier = 0;
        match order {
            0 => {
                b = self.d0(param, first, last);
                local_result = Some(&self.v);
            }
            1 => {
                b = self.d1(param, first, last);
                local_result = Some(&self.dv);
            }
            2 => {
                b = self.d2(param, first, last);
                local_result = Some(&self.d2v);
            }
            _ => {
                b = false;
                local_result = None;
            }
        }
        if !b {
            *ier = order + 1;
        }
        match local_result {
            Some(arr) => {
                // OCCT: for (int ii = 0; ii <= 11; ii++)
                //         AddrResult[ii] = LocalResult[ii];
                for ii in 0..=11usize {
                    result[ii] = arr[ii / 3][ii % 3];
                }
            }
            None => {
                // OCCT copies from the null LocalResult here (undefined
                // behaviour of the default branch, unreachable through the
                // Sweep evaluator which only requests orders 0-2); the
                // Rust form leaves Result untouched once Ier flags it.
            }
        }
    }
}
