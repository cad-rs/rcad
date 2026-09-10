//! OCCT GeomPlate_PlateG0Criterion (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_PlateG0Criterion.cxx / .hxx.
//!
//! GAP leaves: the base class `AdvApp2Var_Criterion` accessors are complete
//! here; the `Value` body's `AdvApp2Var_Patch`/`AdvApp2Var_Context`/
//! `PLib::EvalPoly2Var` machinery still waits for the AdvApp2Var class layer
//! (Node/Iso/Network/Framework/Patch/Context/ApproxAFunc2Var, over the
//! kernel B-spline types).  The full AdvApp2Var ENGINE layer is now landed
//! in [`crate::geomalgo::adv_app2_var`] (SysBase/MathBase subsets incl.
//! mmjacan_/mmjaccv_/mmapcmp_/mmaperx_/mmfmca8_/mmfmca9_/mmfmtb1_/mmmpocur_/
//! mmtrpjj_/mmveps3_/mzsnorm_, ALL ApproxF2var engines mma1*-/mma2can_/
//! mma2cd*-/mma2ce*-/mma2cf*-/mma2er*-/mma2ds*-/mma2fnc_/mmjacpt_ and the
//! MMAPGS*/MLGDRTL/MMJCOBI/MMCMCNP block data); only the class layer above
//! the engines remains.  The criterion carries the base-class members
//! (myMaxValue/myType/myRepartition) and the patch-evaluation entry point
//! preserves the untranslated-dependency failure path; `IsSatisfied` is
//! complete.  The criterion is only consumed through
//! AdvApp2Var_ApproxAFunc2Var, which is itself the MakeApprox GAP.

use glam::{DVec2, DVec3};

/// OCCT AdvApp2Var_CriterionType (AdvApp2Var_CriterionType.hxx L28-32).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdvApp2VarCriterionType {
    /// OCCT AdvApp2Var_Absolute.
    Absolute,
    /// OCCT AdvApp2Var_Relative.
    Relative,
}

/// OCCT AdvApp2Var_CriterionRepartition
/// (AdvApp2Var_CriterionRepartition.hxx L28-32).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdvApp2VarCriterionRepartition {
    /// OCCT AdvApp2Var_Regular.
    Regular,
    /// OCCT AdvApp2Var_Incremental.
    Incremental,
}

/// GAP carrier for OCCT AdvApp2Var_Patch — the patch object whose
/// polynomial-coefficient buffers the criterion evaluates; the
/// AdvApp2Var_Patch class (AdvApp2Var_Patch.cxx, over the not-yet-landed
/// mma2ce*/mma2can_* engines) is the remaining dependency, so only the
/// CritValue slot exists here.
#[derive(Debug, Clone, Copy, Default)]
pub struct AdvApp2VarPatchCarrier {
    crit_value: f64,
}

impl AdvApp2VarPatchCarrier {
    /// OCCT AdvApp2Var_Patch::SetCritValue.
    pub fn set_crit_value(&mut self, v: f64) {
        self.crit_value = v;
    }

    /// OCCT AdvApp2Var_Patch::CritValue.
    pub fn crit_value(&self) -> f64 {
        self.crit_value
    }
}

/// OCCT GeomPlate_PlateG0Criterion (hxx L28-46).
#[derive(Debug, Clone)]
pub struct PlateG0Criterion {
    /// hxx L47-48 (private): NCollection_Sequence<gp_XY> myData.
    my_data: Vec<DVec2>,
    /// hxx L49: NCollection_Sequence<gp_XYZ> myXYZ.
    my_xyz: Vec<DVec3>,
    /// OCCT AdvApp2Var_Criterion base members.
    my_max_value: f64,
    my_type: AdvApp2VarCriterionType,
    my_repartition: AdvApp2VarCriterionRepartition,
}

impl PlateG0Criterion {
    /// OCCT ctor (GeomPlate_PlateG0Criterion.cxx L36-47) — default args
    /// Type = AdvApp2Var_Absolute, Repart = AdvApp2Var_Regular (hxx L38-39).
    pub fn new(
        data: &[DVec2],
        g0data: &[DVec3],
        maximum: f64,
        ctype: AdvApp2VarCriterionType,
        repart: AdvApp2VarCriterionRepartition,
    ) -> Self {
        PlateG0Criterion {
            my_data: data.to_vec(),
            my_xyz: g0data.to_vec(),
            my_max_value: maximum,
            my_type: ctype,
            my_repartition: repart,
        }
    }

    /// OCCT Value(P, C) (GeomPlate_PlateG0Criterion.cxx L51-116).
    ///
    /// GAP leaf: walks the AdvApp2Var_Patch polynomial coefficients through
    /// PLib::EvalPoly2Var against the AdvApp2Var_Context limits — the
    /// AdvApp2Var package is untranslated; the anchor preserves the
    /// dependency failure path.  The body lands together with AdvApp2Var.
    pub fn value(&self, _p: &mut AdvApp2VarPatchCarrier) {
        unimplemented!(
            "GeomPlate_PlateG0Criterion::Value needs AdvApp2Var_Patch/Context + PLib::EvalPoly2Var (untranslated)"
        );
    }

    /// OCCT IsSatisfied(P) (GeomPlate_PlateG0Criterion.cxx L120-123).
    pub fn is_satisfied(&self, p: &AdvApp2VarPatchCarrier) -> bool {
        p.crit_value() < self.my_max_value
    }

    /// OCCT AdvApp2Var_Criterion::MaxValue accessor (base class lxx).
    pub fn max_value(&self) -> f64 {
        self.my_max_value
    }

    /// OCCT AdvApp2Var_Criterion::Type accessor.
    pub fn ctype(&self) -> AdvApp2VarCriterionType {
        self.my_type
    }

    /// OCCT AdvApp2Var_Criterion::Repartition accessor.
    pub fn repartition(&self) -> AdvApp2VarCriterionRepartition {
        self.my_repartition
    }
}
