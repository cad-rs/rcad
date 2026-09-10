//! ChFiDS gap-fill — methods of `ChFiDSSurfData`, `ChFiDS_FaceInterference`
//! and `ChFiDS_CommonPoint` that were missing from `chfi_ds.rs`, translated
//! 1:1 from OCCT TKFillet/ChFiDS.  Kept in a dedicated file because
//! `chfi_ds.rs` exceeds the 2000-line guideline.

use glam::{DVec2, DVec3};
use rcad_kernel::topo::topods::Orientation;

use super::chfi_ds::{
    ChFiDS_CommonPoint, ChFiDS_FaceInterference, ChFiDSCircSectionArray, ChFiDSSurfData,
};

// =========================================================================
// OCCT ChFiDS_SurfData — missing method translations
// =========================================================================

impl ChFiDSSurfData {
    /// OCCT ChFiDS_SurfData.cxx L44-72 — Copy(Other).  OCCT also assigns
    /// `simul = Other->simul`; the simul slot is a pending boundary in rcad
    /// (handle<Standard_Transient> has no rcad equivalent yet).
    pub fn copy(&mut self, other: &ChFiDSSurfData) {
        self.index_of_s1 = other.index_of_s1;
        self.index_of_s2 = other.index_of_s2;
        self.index_of_conge = other.index_of_conge;
        self.orientation = other.orientation;
        self.intf1 = other.intf1.clone();
        self.intf2 = other.intf2.clone();

        self.pfirst_on_s1 = other.pfirst_on_s1.clone();
        self.plast_on_s1 = other.plast_on_s1.clone();
        self.pfirst_on_s2 = other.pfirst_on_s2.clone();
        self.plast_on_s2 = other.plast_on_s2.clone();

        self.ufspine = other.ufspine;
        self.ulspine = other.ulspine;

        // OCCT: simul = Other->simul.
        self.simul = other.simul.clone();

        self.p2df1 = other.p2df1;
        self.p2dl1 = other.p2dl1;
        self.p2df2 = other.p2df2;
        self.p2dl2 = other.p2dl2;

        self.myfirstextend = other.myfirstextend;
        self.mylastextend = other.mylastextend;

        self.twistons1 = other.twistons1;
        self.twistons2 = other.twistons2;
    }

    /// OCCT ChFiDS_SurfData.lxx L19-22 — IndexOfS1().
    pub fn index_of_s1(&self) -> i32 {
        self.index_of_s1
    }

    /// OCCT ChFiDS_SurfData.lxx L50-53 — IndexOfS2().
    pub fn index_of_s2(&self) -> i32 {
        self.index_of_s2
    }

    /// OCCT ChFiDS_SurfData.lxx L33-38 — IndexOfC1().
    pub fn index_of_c1(&self) -> i32 {
        if !self.isoncurv1 {
            panic!("Standard_Failure: Interference pas sur courbe");
        }
        self.index_of_c1
    }

    /// OCCT ChFiDS_SurfData.lxx L42-46 — SetIndexOfC1(theIndex).
    pub fn set_index_of_c1(&mut self, the_index: i32) {
        self.index_of_c1 = the_index;
        self.isoncurv1 = the_index != 0;
    }

    /// OCCT ChFiDS_SurfData.lxx L64-69 — IndexOfC2().
    pub fn index_of_c2(&self) -> i32 {
        if !self.isoncurv2 {
            panic!("Standard_Failure: Interference pas sur courbe");
        }
        self.index_of_c2
    }

    /// OCCT ChFiDS_SurfData.lxx L73-77 — SetIndexOfC2(theIndex).
    pub fn set_index_of_c2(&mut self, the_index: i32) {
        self.index_of_c2 = the_index;
        self.isoncurv2 = the_index != 0;
    }

    /// OCCT ChFiDS_SurfData.lxx L179-182 — ChangeVertexFirstOnS1().
    pub fn change_vertex_first_on_s1(&mut self) -> &mut ChFiDS_CommonPoint {
        &mut self.pfirst_on_s1
    }

    /// OCCT ChFiDS_SurfData.lxx L186-189 — ChangeVertexLastOnS1().
    pub fn change_vertex_last_on_s1(&mut self) -> &mut ChFiDS_CommonPoint {
        &mut self.plast_on_s1
    }

    /// OCCT ChFiDS_SurfData.lxx L193-196 — ChangeVertexFirstOnS2().
    pub fn change_vertex_first_on_s2(&mut self) -> &mut ChFiDS_CommonPoint {
        &mut self.pfirst_on_s2
    }

    /// OCCT ChFiDS_SurfData.lxx L200-203 — ChangeVertexLastOnS2().
    pub fn change_vertex_last_on_s2(&mut self) -> &mut ChFiDS_CommonPoint {
        &mut self.plast_on_s2
    }

    /// OCCT ChFiDS_SurfData.lxx L207-212 — IsOnCurve(OnS).
    pub fn is_on_curve(&self, on_s: i32) -> bool {
        if on_s == 1 {
            self.isoncurv1
        } else {
            self.isoncurv2
        }
    }

    /// OCCT ChFiDS_SurfData.lxx L216-227 — IndexOfC(OnS).
    pub fn index_of_c(&self, on_s: i32) -> i32 {
        if on_s == 1 {
            if !self.isoncurv1 {
                panic!("Standard_Failure: Interference pas sur courbe");
            }
            return self.index_of_c1;
        }
        if !self.isoncurv2 {
            panic!("Standard_Failure: Interference pas sur courbe");
        }
        self.index_of_c2
    }

    /// OCCT ChFiDS_SurfData.lxx L231-234 — TwistOnS1().
    pub fn twist_on_s1(&self) -> bool {
        self.twistons1
    }

    /// OCCT ChFiDS_SurfData.lxx L236-239 — TwistOnS2().
    pub fn twist_on_s2(&self) -> bool {
        self.twistons2
    }

    /// OCCT ChFiDS_SurfData.lxx L241-244 — TwistOnS1(T).
    pub fn set_twist_on_s1(&mut self, t: bool) {
        self.twistons1 = t;
    }

    /// OCCT ChFiDS_SurfData.lxx L246-249 — TwistOnS2(T).
    pub fn set_twist_on_s2(&mut self, t: bool) {
        self.twistons2 = t;
    }

    /// OCCT ChFiDS_SurfData.cxx L190-193 — FirstExtensionValue().
    pub fn first_extension_value(&self) -> f64 {
        self.myfirstextend
    }

    /// OCCT ChFiDS_SurfData.cxx L197-200 — LastExtensionValue().
    pub fn last_extension_value(&self) -> f64 {
        self.mylastextend
    }

    /// OCCT ChFiDS_SurfData.cxx L204-207 — FirstExtensionValue(Extend).
    pub fn set_first_extension_value(&mut self, extend: f64) {
        self.myfirstextend = extend;
    }

    /// OCCT ChFiDS_SurfData.cxx L211-214 — LastExtensionValue(Extend).
    pub fn set_last_extension_value(&mut self, extend: f64) {
        self.mylastextend = extend;
    }

    /// OCCT ChFiDS_SurfData.cxx L239-248 — Get2dPoints(P2df1..P2dl2).
    pub fn get_2d_points(&self) -> (DVec2, DVec2, DVec2, DVec2) {
        (self.p2df1, self.p2dl1, self.p2df2, self.p2dl2)
    }

    /// OCCT ChFiDS_SurfData.cxx L252-268 — Get2dPoints(First, OnS).
    pub fn get_2d_point(&self, first: bool, on_s: i32) -> DVec2 {
        if first && on_s == 1 {
            self.p2df1
        } else if !first && on_s == 1 {
            self.p2dl1
        } else if first && on_s == 2 {
            self.p2df2
        } else {
            self.p2dl2
        }
    }

    /// OCCT ChFiDS_SurfData.cxx L272-281 — Set2dPoints(P2df1..P2dl2).
    pub fn set_2d_points(&mut self, p2df1: DVec2, p2dl1: DVec2, p2df2: DVec2, p2dl2: DVec2) {
        self.p2df1 = p2df1;
        self.p2dl1 = p2dl1;
        self.p2df2 = p2df2;
        self.p2dl2 = p2dl2;
    }

    /// OCCT ChFiDS_SurfData.cxx L218-221 — Simul().  Architecture mapping:
    /// the transient-handle down-cast to NCollection_HArray1<ChFiDS_CircSection>
    /// (FilBuilder Sect, ChFi3d_FilBuilder.cxx L465) is implicit — the rcad
    /// slot already carries the concrete array.
    pub fn simul(&self) -> Option<ChFiDSCircSectionArray> {
        self.simul.clone()
    }

    /// OCCT ChFiDS_SurfData.cxx L225-228 — SetSimul(S).  The handle
    /// parameter may be null (SimulKPart calls SetSimul unconditionally,
    /// ChFi3d_FilBuilder.cxx L559), so the rcad slot takes the Option.
    pub fn set_simul(&mut self, s: Option<ChFiDSCircSectionArray>) {
        self.simul = s;
    }

    /// OCCT ChFiDS_SurfData.cxx L232-235 — ResetSimul().
    pub fn reset_simul(&mut self) {
        self.simul = None;
    }
}

// =========================================================================
// OCCT ChFiDS_FaceInterference — missing method translations
// =========================================================================

impl ChFiDS_FaceInterference {
    /// OCCT ChFiDS_FaceInterference.lxx L32-35 — SetLineIndex(I).
    pub fn set_line_index(&mut self, i: i32) {
        self.lineindex = i;
    }

    /// OCCT ChFiDS_FaceInterference.lxx L39-42 — SetFirstParameter(U1).
    pub fn set_first_parameter(&mut self, u1: f64) {
        self.first_param = u1;
    }

    /// OCCT ChFiDS_FaceInterference.lxx L46-49 — SetLastParameter(U1).
    pub fn set_last_parameter(&mut self, u1: f64) {
        self.last_param = u1;
    }

    /// OCCT ChFiDS_FaceInterference.cxx L42-45 — SetTransition(Trans).
    pub fn set_transition(&mut self, trans: Orientation) {
        self.line_transition = trans;
    }
}

// =========================================================================
// OCCT ChFiDS_CommonPoint — missing method translations
// =========================================================================

impl ChFiDS_CommonPoint {
    /// OCCT ChFiDS_CommonPoint.hxx L127 — HasVector().
    pub fn has_vector(&self) -> bool {
        self.hasvector
    }
}

#[cfg(test)]
mod simul_tests {
    use super::super::chfi_ds::{ChFiDSCircSection, ChFiDSSurfData};

    /// OCCT anchor: Simul / SetSimul / ResetSimul (ChFiDS_SurfData.cxx
    /// L218-235) and the simul transfer in Copy (L44-72).
    #[test]
    fn simul_slot_roundtrip() {
        let mut sd = ChFiDSSurfData::default();
        assert!(sd.simul().is_none());

        let sec = vec![ChFiDSCircSection::new(), ChFiDSCircSection::new()];
        sd.set_simul(Some(sec));
        let got = sd.simul().expect("simul stored");
        assert_eq!(got.len(), 2);

        let mut copied = ChFiDSSurfData::default();
        copied.copy(&sd);
        assert!(copied.simul().is_some(), "Copy transfers the simul slot");

        sd.reset_simul();
        assert!(sd.simul().is_none());
    }
}
