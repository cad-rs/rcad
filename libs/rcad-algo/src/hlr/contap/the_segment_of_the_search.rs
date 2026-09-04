// OCCT Contap_TheSegmentOfTheSearch (TKHLR) — a part of a domain arc
// solution of Contap_TheSearch.
//
// Contap_TheSegmentOfTheSearch.hxx L27-112 + _0.cxx L24-43.

use super::point::Arc;
use super::the_path_point_of_the_search::ThePathPointOfTheSearch;

/// OCCT Contap_TheSegmentOfTheSearch.
#[derive(Clone)]
pub struct TheSegmentOfTheSearch {
    arc: Option<Arc>,
    hasfp: bool,
    thefp: ThePathPointOfTheSearch,
    haslp: bool,
    thelp: ThePathPointOfTheSearch,
}

impl TheSegmentOfTheSearch {
    /// OCCT Contap_TheSegmentOfTheSearch() (_0.cxx L24-28).
    pub fn new() -> Self {
        TheSegmentOfTheSearch {
            arc: None,
            hasfp: false,
            thefp: ThePathPointOfTheSearch::new(),
            haslp: false,
            thelp: ThePathPointOfTheSearch::new(),
        }
    }

    /// OCCT SetValue(A) (hxx L72-77) — defines the concerned arc.
    pub fn set_value(&mut self, a: Arc) {
        self.hasfp = false;
        self.haslp = false;
        self.arc = Some(a);
    }

    /// OCCT SetLimitPoint(V, First) (_0.cxx L30-43).
    pub fn set_limit_point(&mut self, v: &ThePathPointOfTheSearch, first: bool) {
        if first {
            self.hasfp = true;
            self.thefp = v.clone();
        } else {
            self.haslp = true;
            self.thelp = v.clone();
        }
    }

    /// OCCT Curve (hxx L79-82).
    pub fn curve(&self) -> Arc {
        self.arc.clone().expect("TheSegmentOfTheSearch::Curve")
    }

    /// OCCT HasFirstPoint (hxx L84-87).
    pub fn has_first_point(&self) -> bool {
        self.hasfp
    }

    /// OCCT FirstPoint (hxx L89-96).
    pub fn first_point(&self) -> &ThePathPointOfTheSearch {
        if !self.hasfp {
            panic!("Standard_DomainError: TheSegmentOfTheSearch::FirstPoint");
        }
        &self.thefp
    }

    /// OCCT HasLastPoint (hxx L98-101).
    pub fn has_last_point(&self) -> bool {
        self.haslp
    }

    /// OCCT LastPoint (hxx L103-110).
    pub fn last_point(&self) -> &ThePathPointOfTheSearch {
        if !self.haslp {
            panic!("Standard_DomainError: TheSegmentOfTheSearch::LastPoint");
        }
        &self.thelp
    }
}

impl Default for TheSegmentOfTheSearch {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::{DVec2, DVec3};
    use rcad_kernel::geom::Curve2d;

    fn test_arc() -> Arc {
        std::sync::Arc::new(Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::new(1.0, 0.0),
        }))
    }

    /// OCCT anchor: SetLimitPoint stores first/last independently;
    /// FirstPoint/LastPoint raise DomainError when absent (hxx L89-110).
    #[test]
    fn segment_limit_points() {
        let mut seg = TheSegmentOfTheSearch::new();
        seg.set_value(test_arc());
        assert!(!seg.has_first_point());
        assert!(!seg.has_last_point());

        let p = ThePathPointOfTheSearch::without_vertex(DVec3::ONE, 1e-5, test_arc(), 0.25);
        seg.set_limit_point(&p, true);
        assert!(seg.has_first_point());
        assert!((seg.first_point().parameter() - 0.25).abs() < 1e-15);
        assert!(!seg.has_last_point());
    }
}
