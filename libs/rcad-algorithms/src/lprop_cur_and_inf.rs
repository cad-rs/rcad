//! LProp-style curve curvature extremum and inflection point storage.
//!
//! LProp_CurAndInf — stores curvature extrema (min/max)
//!                  and inflection points of a curve, sorted by parameter.

/// Type of curvature/inflection point (OCCT LProp_CIType).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CIType {
    Inflection,
    MinCur,
    MaxCur,
}

/// Storage for curvature extrema and inflection points.
///
/// LProp_CurAndInf. Points are stored sorted by parameter.
#[derive(Debug, Clone)]
pub struct CurAndInf {
    params: Vec<f64>,
    types: Vec<CIType>,
}

impl CurAndInf {
    pub fn new() -> Self { Self { params: Vec::new(), types: Vec::new() } }

    pub fn is_empty(&self) -> bool { self.params.is_empty() }
    pub fn nb_points(&self) -> usize { self.params.len() }

    pub fn add_inflection(&mut self, param: f64) {
        self.params.push(param);
        self.types.push(CIType::Inflection);
        self.sort();
    }

    pub fn add_ext_cur(&mut self, param: f64, is_min: bool) {
        self.params.push(param);
        self.types.push(if is_min { CIType::MinCur } else { CIType::MaxCur });
        self.sort();
    }

    pub fn clear(&mut self) { self.params.clear(); self.types.clear(); }

    pub fn parameter(&self, idx: usize) -> f64 {
        self.params[idx]
    }

    pub fn ci_type(&self, idx: usize) -> CIType {
        self.types[idx]
    }

    fn sort(&mut self) {
        let mut indices: Vec<usize> = (0..self.params.len()).collect();
        indices.sort_by(|&a, &b| self.params[a].partial_cmp(&self.params[b]).unwrap());
        self.params = indices.iter().map(|&i| self.params[i]).collect();
        self.types = indices.iter().map(|&i| self.types[i]).collect();
    }
}

// =============================================================================
// Tests — translated from LProp_CurAndInf_Test.cxx
// =============================================================================


