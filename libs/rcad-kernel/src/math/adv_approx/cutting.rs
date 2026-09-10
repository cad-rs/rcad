//! OCCT AdvApprox_Cutting + AdvApprox_DichoCutting (TKG3d).

/// OCCT AdvApprox_Cutting — method to split an interval in two.  Returns
/// `Some(cuttingvalue)` for "Large" (cutting allowed) and `None` otherwise;
/// the value is computed in both cases by the OCCT implementations but only
/// consumed when the cut happens, so the collapsed `Option` is faithful.
pub trait Cutting {
    /// OCCT `Value(a, b, cuttingvalue)`.
    fn value(&self, a: f64, b: f64) -> Option<f64>;
}

/// OCCT AdvApprox_DichoCutting — cuts at (a + b) / 2.
pub struct DichoCutting;

impl Cutting for DichoCutting {
    fn value(&self, a: f64, b: f64) -> Option<f64> {
        // Minimum length of an interval for F(U,V): EPS1=1.e-9 (cf. MEPS1).
        let lgmin = 10.0 * crate::core::precision::p_confusion();
        let cuttingvalue = (a + b) / 2.0;
        if (b - a).abs() >= 2.0 * lgmin {
            Some(cuttingvalue)
        } else {
            None
        }
    }
}
