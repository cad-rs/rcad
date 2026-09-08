//! OCCT MAT_Edge — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_Edge.hxx L17-65, MAT_Edge.cxx L23-80

use std::sync::{Arc, RwLock};

use super::mat_bisector::HandleMatBisector;

/// OCCT `occ::handle<MAT_Edge>`
pub type HandleMatEdge = Arc<RwLock<MatEdge>>;

/// OCCT MAT_Edge.hxx L27-63
pub struct MatEdge {
    theedgenumber: i32,
    thefirstbisector: Option<HandleMatBisector>,
    thesecondbisector: Option<HandleMatBisector>,
    thedistance: f64,
    theintersectionpoint: i32,
}

impl MatEdge {
    /// OCCT MAT_Edge.cxx L23-28
    pub fn new() -> Self {
        MatEdge {
            theedgenumber: 0,
            thefirstbisector: None,
            thesecondbisector: None,
            thedistance: 0.0,
            theintersectionpoint: 0,
        }
    }

    /// OCCT MAT_Edge.cxx L30-33 (setter overload of EdgeNumber)
    pub fn set_edge_number(&mut self, anumber: i32) {
        self.theedgenumber = anumber;
    }

    /// OCCT MAT_Edge.cxx L35-38
    pub fn set_first_bisector(&mut self, abisector: &HandleMatBisector) {
        self.thefirstbisector = Some(abisector.clone());
    }

    /// OCCT MAT_Edge.cxx L40-43
    pub fn set_second_bisector(&mut self, abisector: &HandleMatBisector) {
        self.thesecondbisector = Some(abisector.clone());
    }

    // OCCT ~MAT2d_Mat2d (cxx L1878-1910): the destructor nulls the bisector
    // handles — the clear entry points the MAT2d Drop path calls.
    pub fn clear_first_bisector(&mut self) {
        self.thefirstbisector = None;
    }

    pub fn clear_second_bisector(&mut self) {
        self.thesecondbisector = None;
    }

    /// OCCT MAT_Edge.cxx L45-48
    pub fn set_distance(&mut self, adistance: f64) {
        self.thedistance = adistance;
    }

    /// OCCT MAT_Edge.cxx L50-53
    pub fn set_intersection_point(&mut self, apoint: i32) {
        self.theintersectionpoint = apoint;
    }

    /// OCCT MAT_Edge.cxx L55-58 (getter overload of EdgeNumber)
    pub fn edge_number(&self) -> i32 {
        self.theedgenumber
    }

    /// OCCT MAT_Edge.cxx L60-63
    pub fn first_bisector(&self) -> Option<HandleMatBisector> {
        self.thefirstbisector.clone()
    }

    /// OCCT MAT_Edge.cxx L65-68
    pub fn second_bisector(&self) -> Option<HandleMatBisector> {
        self.thesecondbisector.clone()
    }

    /// OCCT MAT_Edge.cxx L70-73
    pub fn distance(&self) -> f64 {
        self.thedistance
    }

    /// OCCT MAT_Edge.cxx L75-78
    pub fn intersection_point(&self) -> i32 {
        self.theintersectionpoint
    }

    /// OCCT MAT_Edge.cxx L80
    pub fn dump(&self, _ashift: i32, _alevel: i32) {}
}

impl Default for MatEdge {
    fn default() -> Self {
        Self::new()
    }
}
