//! OCCT BOPAlgo_SectionAttribute (BOPAlgo_SectionAttribute.hxx).
//! Section attributes for boolean operations (method, tolerance, etc.).
//!
//! Controls whether pcurves are stored on section edges for each face.

/// OCCT BOPAlgo_SectionAttribute.
#[derive(Debug, Clone)]
pub struct SectionAttribute {
    pub method: i32,
    pub section_tolerance: f64,
    // OCCT BOPAlgo_SectionAttribute.hxx L24-38: myApproximation, myPCurve1, myPCurve2
    // all default to true.
    pub approximation: bool, // myApproximation
    pub pcurve_on_s1: bool,  // myPCurve1
    pub pcurve_on_s2: bool,  // myPCurve2
}

impl Default for SectionAttribute {
    fn default() -> Self {
        Self {
            method: 0,
            section_tolerance: 1e-7,
            approximation: true,
            pcurve_on_s1: true,
            pcurve_on_s2: true,
        }
    }
}
