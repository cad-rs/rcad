//! LProp-style curve curvature extremum and inflection point storage.
//!
//! ✅ OCCT-aligned: LProp_CurAndInf — stores curvature extrema (min/max)
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
/// OCCT-aligned: LProp_CurAndInf. Points are stored sorted by parameter.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_constructor_is_empty() {
        let r = CurAndInf::new();
        assert!(r.is_empty());
        assert_eq!(r.nb_points(), 0);
    }

    #[test]
    fn add_inflection() {
        let mut r = CurAndInf::new();
        r.add_inflection(1.5);
        assert!(!r.is_empty());
        assert_eq!(r.nb_points(), 1);
        assert!((r.parameter(0) - 1.5).abs() < 1e-12);
        assert_eq!(r.ci_type(0), CIType::Inflection);
    }

    #[test]
    fn add_ext_cur_minimum() {
        let mut r = CurAndInf::new();
        r.add_ext_cur(2.0, true);
        assert_eq!(r.nb_points(), 1);
        assert_eq!(r.ci_type(0), CIType::MinCur);
    }

    #[test]
    fn add_ext_cur_maximum() {
        let mut r = CurAndInf::new();
        r.add_ext_cur(3.0, false);
        assert_eq!(r.nb_points(), 1);
        assert_eq!(r.ci_type(0), CIType::MaxCur);
    }

    #[test]
    fn multiple_points_sorted_by_parameter() {
        let mut r = CurAndInf::new();
        r.add_inflection(5.0);
        r.add_ext_cur(1.0, true);
        r.add_ext_cur(3.0, false);
        r.add_inflection(7.0);
        assert_eq!(r.nb_points(), 4);
        for i in 1..r.nb_points() {
            assert!(r.parameter(i - 1) <= r.parameter(i));
        }
    }

    #[test]
    fn clear_and_refill() {
        let mut r = CurAndInf::new();
        r.add_inflection(1.0);
        r.clear();
        r.add_ext_cur(5.0, false);
        assert_eq!(r.nb_points(), 1);
        assert_eq!(r.ci_type(0), CIType::MaxCur);
    }
}
