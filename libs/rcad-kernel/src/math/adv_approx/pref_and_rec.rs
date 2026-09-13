//! OCCT AdvApprox_PrefAndRec (TKG3d/AdvApprox) — the cutting tool carrying a
//! list of preferential points (pi)i and a list of recommended points.
//!
//! 1:1 translation of `AdvApprox_PrefAndRec.hxx` (L28-59) and
//! `AdvApprox_PrefAndRec.cxx` (L22-75).

use super::cutting::Cutting;
use crate::core::precision;

/// OCCT AdvApprox_PrefAndRec (hxx L33-59) — inherits class Cutting; contains
/// a list of preferential points (pi)i and a list of Recommended points used
/// in cutting management.  If Cutting is necessary in [a,b], we cut at the di
/// nearest from (a+b)/2.
pub struct PrefAndRec {
    /// OCCT: NCollection_Array1<double> myRecCutting.
    my_rec_cutting: Vec<f64>,
    /// OCCT: NCollection_Array1<double> myPrefCutting.
    my_pref_cutting: Vec<f64>,
    /// OCCT: double myWeight.
    my_weight: f64,
}

impl PrefAndRec {
    /// OCCT AdvApprox_PrefAndRec(RecCut, PrefCut, Weight=5)
    /// (cxx L32-46): the arrays are copied and the weight must be > 1.
    pub fn new(rec_cut: &[f64], pref_cut: &[f64], weight: f64) -> Self {
        let this = PrefAndRec {
            my_rec_cutting: rec_cut.to_vec(),
            my_pref_cutting: pref_cut.to_vec(),
            my_weight: weight,
        };
        if this.my_weight <= 1.0 {
            // OCCT: throw Standard_DomainError("PrefAndRec : Weight is too small").
            panic!("Standard_DomainError: PrefAndRec : Weight is too small");
        }
        this
    }

    /// OCCT AdvApprox_PrefAndRec(RecCut, PrefCut) — the default Weight = 5.
    pub fn with_default_weight(rec_cut: &[f64], pref_cut: &[f64]) -> Self {
        Self::new(rec_cut, pref_cut, 5.0)
    }
}

impl Default for PrefAndRec {
    fn default() -> Self {
        Self::with_default_weight(&[], &[])
    }
}

impl Cutting for PrefAndRec {
    /// OCCT AdvApprox_PrefAndRec::Value(a, b, cuttingvalue) (cxx L48-75).
    ///
    /// cutting value is
    /// - the recommended point nearest of (a+b)/2
    ///   if pi is in ]a,b[ or else
    /// - the preferential point nearest of (a+b)/2
    ///   if pi is in ](r*a+b)/(r+1), (a+r*b)/(r+1)[ where r = Weight
    /// - or (a+b)/2 else.
    ///
    /// The bool result ("Large": the cut is at least lgmin from both bounds)
    /// is the `Option` discriminant of the `Cutting` trait — the OCCT
    /// cuttingvalue written on the `false` return is never consumed by
    /// `Approximation` (AdvApprox_ApproxAFunction.cxx L543-547).
    fn value(&self, a: f64, b: f64) -> Option<f64> {
        // Minimum length of a parametric interval: 10*PConfusion() (cxx L50).
        let lgmin = 10.0 * precision::p_confusion();
        let mil = (a + b) / 2.0;
        let mut dist;
        let mut is_found = false;

        let mut cut = mil;

        // Search for a preferred cutting point (cxx L56-66).
        dist = ((a * self.my_weight + b) / (1.0 + self.my_weight) - mil).abs();
        for i in 1..=self.my_pref_cutting.len() {
            let pi = self.my_pref_cutting[i - 1];
            if dist > (mil - pi).abs() {
                cut = pi;
                dist = (mil - cut).abs();
                is_found = true;
            }
        }

        // Search for a recommended cutting point (cxx L68-78).
        if !is_found {
            dist = ((a - b) / 2.0).abs();
            for i in 1..=self.my_rec_cutting.len() {
                let ri = self.my_rec_cutting[i - 1];
                if (dist - lgmin) > (mil - ri).abs() {
                    cut = ri;
                    dist = (mil - cut).abs();
                }
            }
        }

        // Result (cxx L80-82).
        if (cut - a).abs() >= lgmin && (b - cut).abs() >= lgmin {
            Some(cut)
        } else {
            None
        }
    }
}
