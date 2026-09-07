//! ChFiDS gap-fill — methods of `ChFiDS_Stripe`, `ChFiDSStripeMap` and
//! `ChFiDSMap` that were missing from `chfi_ds.rs`, translated 1:1 from
//! OCCT TKFillet/ChFiDS.  Kept in a dedicated file because `chfi_ds.rs`
//! exceeds the 2000-line guideline.

use rcad_kernel::topo::topods::{Orientation, Shape};

use super::chfi_ds::{ChFiDSMap, ChFiDSStripe, ChFiDSStripeMap, SharedStripe, SharedSurfData};

// =========================================================================
// OCCT ChFiDS_Stripe — missing method translations
// =========================================================================

impl ChFiDSStripe {
    /// OCCT ChFiDS_Stripe.cxx L57-69 — Parameters(First, Pdeb, Pfin).
    pub fn parameters(&self, first: bool) -> (f64, f64) {
        if first {
            (self.pardeb1, self.parfin1)
        } else {
            (self.pardeb2, self.parfin2)
        }
    }

    /// OCCT ChFiDS_Stripe.lxx L107-111 — ChangeFirstParameters(Pdeb, Pfin).
    pub fn change_first_parameters(&mut self, pdeb: f64, pfin: f64) {
        self.pardeb1 = pdeb;
        self.parfin1 = pfin;
    }

    /// OCCT ChFiDS_Stripe.lxx L115-119 — ChangeLastParameters(Pdeb, Pfin).
    pub fn change_last_parameters(&mut self, pdeb: f64, pfin: f64) {
        self.pardeb2 = pdeb;
        self.parfin2 = pfin;
    }

    /// OCCT ChFiDS_Stripe.cxx L89-99 — Curve(First).
    pub fn curve(&self, first: bool) -> i32 {
        if first {
            self.index_ofcurve1
        } else {
            self.index_ofcurve2
        }
    }

    /// OCCT ChFiDS_Stripe.lxx L137-140 — ChangeFirstCurve(Index).
    pub fn change_first_curve(&mut self, index: i32) {
        self.index_ofcurve1 = index;
    }

    /// OCCT ChFiDS_Stripe.lxx L144-147 — ChangeLastCurve(Index).
    pub fn change_last_curve(&mut self, index: i32) {
        self.index_ofcurve2 = index;
    }

    /// OCCT ChFiDS_Stripe.lxx L19-23 — SetOfSurfData().
    pub fn set_of_surf_data(&self) -> &Vec<SharedSurfData> {
        &self.my_hdata
    }

    /// OCCT ChFiDS_Stripe.lxx L55-59 — ChangeSetOfSurfData().
    pub fn change_set_of_surf_data(&mut self) -> &mut Vec<SharedSurfData> {
        &mut self.my_hdata
    }

    /// OCCT ChFiDS_Stripe.lxx L34-37 — OrientationOnFace1().
    pub fn orientation_on_face1(&self) -> Orientation {
        self.my_or1
    }

    /// OCCT ChFiDS_Stripe.lxx L41-44 — OrientationOnFace2().
    pub fn orientation_on_face2(&self) -> Orientation {
        self.my_or2
    }

    /// OCCT ChFiDS_Stripe.lxx L70-73 — OrientationOnFace1(Or1).
    pub fn set_orientation_on_face1(&mut self, or1: Orientation) {
        self.my_or1 = or1;
    }

    /// OCCT ChFiDS_Stripe.lxx L77-80 — OrientationOnFace2(Or2).
    pub fn set_orientation_on_face2(&mut self, or2: Orientation) {
        self.my_or2 = or2;
    }

    /// OCCT ChFiDS_Stripe.lxx L48-51 — Choix().
    pub fn choix(&self) -> i32 {
        self.my_choix
    }

    /// OCCT ChFiDS_Stripe.lxx L84-87 — Choix(C).
    pub fn set_choix(&mut self, c: i32) {
        self.my_choix = c;
    }

    /// OCCT ChFiDS_Stripe.lxx L179-182 — IndexFirstPointOnS1().
    pub fn index_first_point_on_s1(&self) -> i32 {
        self.indexfirst_pon_s1
    }

    /// OCCT ChFiDS_Stripe.lxx L186-189 — IndexLastPointOnS1().
    pub fn index_last_point_on_s1(&self) -> i32 {
        self.indexlast_pon_s1
    }

    /// OCCT ChFiDS_Stripe.lxx L193-196 — IndexFirstPointOnS2().
    pub fn index_first_point_on_s2(&self) -> i32 {
        self.indexfirst_pon_s2
    }

    /// OCCT ChFiDS_Stripe.lxx L200-203 — IndexLastPointOnS2().
    pub fn index_last_point_on_s2(&self) -> i32 {
        self.indexlast_pon_s2
    }

    /// OCCT ChFiDS_Stripe.lxx L207-210 — ChangeIndexFirstPointOnS1(Index).
    pub fn change_index_first_point_on_s1(&mut self, index: i32) {
        self.indexfirst_pon_s1 = index;
    }

    /// OCCT ChFiDS_Stripe.lxx L214-217 — ChangeIndexLastPointOnS1(Index).
    pub fn change_index_last_point_on_s1(&mut self, index: i32) {
        self.indexlast_pon_s1 = index;
    }

    /// OCCT ChFiDS_Stripe.lxx L221-224 — ChangeIndexFirstPointOnS2(Index).
    pub fn change_index_first_point_on_s2(&mut self, index: i32) {
        self.indexfirst_pon_s2 = index;
    }

    /// OCCT ChFiDS_Stripe.lxx L228-231 — ChangeIndexLastPointOnS2(Index).
    pub fn change_index_last_point_on_s2(&mut self, index: i32) {
        self.indexlast_pon_s2 = index;
    }

    /// OCCT ChFiDS_Stripe.lxx L249-252 — FirstPCurveOrientation(O).
    pub fn set_first_pcurve_orientation(&mut self, o: Orientation) {
        self.orcurv1 = o;
    }

    /// OCCT ChFiDS_Stripe.lxx L256-259 — LastPCurveOrientation(O).
    pub fn set_last_pcurve_orientation(&mut self, o: Orientation) {
        self.orcurv2 = o;
    }
}

// =========================================================================
// OCCT ChFiDS_StripeMap — missing method translations
// =========================================================================

impl ChFiDSStripeMap {
    /// OCCT ChFiDS_StripeMap.cxx L39-43 — FindFromKey(V).
    pub fn find_from_key(&self, v: &Shape) -> &Vec<SharedStripe> {
        self.my_map.get(&v.ptr_id()).expect("stripe map key")
    }
}

// =========================================================================
// OCCT ChFiDS_Map — missing method translations
// =========================================================================

impl ChFiDSMap {
    /// OCCT ChFiDS_Map.cxx L48-51 — FindFromIndex(I) (1-based).
    pub fn find_from_index(&self, i: usize) -> &Vec<Shape> {
        self.my_map
            .get(&self.my_keys[i - 1].ptr_id())
            .expect("map index")
    }
}
