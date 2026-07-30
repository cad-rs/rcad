//! OCCT BOPAlgo_SectionAttribute (BOPAlgo_SectionAttribute.hxx).
//! Section attributes for boolean operations (method, tolerance, etc.).
//!
//! Controls whether pcurves are stored on section edges for each face.

/// OCCT BOPAlgo_SectionAttribute.
#[derive(Debug, Clone)]
pub struct SectionAttribute {
    pub method: i32,
    pub section_tolerance: f64,
}

impl Default for SectionAttribute {
    fn default() -> Self {
        Self { method: 0, section_tolerance: 1e-7 }
    }
}
