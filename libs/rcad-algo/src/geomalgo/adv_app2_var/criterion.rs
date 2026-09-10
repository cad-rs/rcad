//! OCCT AdvApp2Var_Criterion (AdvApp2Var_Criterion.hxx + .cxx),
//! AdvApp2Var_CriterionType.hxx and AdvApp2Var_CriterionRepartition.hxx.
//!
//! Encoding note (architecture): the OCCT abstract base class with protected
//! data members (myMaxValue / myType / myRepartition) becomes a trait; each
//! concrete criterion (GeomPlate_PlateG0Criterion, ...) keeps the three base
//! members and serves the accessor methods.  The virtual destructor and the
//! DEFINE_STANDARD_ALLOC machinery have no rcad counterpart.

/// OCCT AdvApp2Var_CriterionType (AdvApp2Var_CriterionType.hxx L25-29) -
/// influence of the criterion on cutting process.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriterionType {
    /// OCCT AdvApp2Var_Absolute - cutting when criterion is not satisfied.
    Absolute,
    /// OCCT AdvApp2Var_Relative - deactivation of the compute of the error
    /// max / cutting when error max is not good or if error max is good and
    /// criterion is not satisfied.
    Relative,
}

/// OCCT AdvApp2Var_CriterionRepartition
/// (AdvApp2Var_CriterionRepartition.hxx L25-29) - way of cutting process.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriterionRepartition {
    /// OCCT AdvApp2Var_Regular - all new cutting points at each step of
    /// cutting process.
    Regular,
    /// OCCT AdvApp2Var_Incremental - add one new cutting point at each step
    /// of cutting process.
    Incremental,
}

/// OCCT AdvApp2Var_Criterion (AdvApp2Var_Criterion.hxx L32-53) - this class
/// contains a given criterion to be satisfied.
pub trait Criterion {
    /// OCCT Value(P, C) (pure virtual, hxx L39).
    fn value(&self, p: &mut super::patch::Patch, c: &super::context::Context);

    /// OCCT IsSatisfied(P) (pure virtual, hxx L41).
    fn is_satisfied(&self, p: &super::patch::Patch) -> bool;

    /// OCCT MaxValue() (AdvApp2Var_Criterion.cxx L23-26) - reads the
    /// protected member myMaxValue (hxx L50), held by the implementer.
    fn max_value(&self) -> f64;

    /// OCCT Type() (AdvApp2Var_Criterion.cxx L30-33) - reads the protected
    /// member myType (hxx L51), held by the implementer.
    fn crit_type(&self) -> CriterionType;

    /// OCCT Repartition() (AdvApp2Var_Criterion.cxx L37-40) - reads the
    /// protected member myRepartition (hxx L52), held by the implementer.
    fn repartition(&self) -> CriterionRepartition;
}
