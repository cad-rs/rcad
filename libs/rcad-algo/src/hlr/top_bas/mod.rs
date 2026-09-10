//! OCCT TopBas package (TKHLR) — basic interference definitions.
//!
//! 1:1 translation of `TopBas_TestInterference.hxx` (L27-158, fully inline)
//! + `TopBas_TestInterference_0.cxx` (L17-36).

use rcad_kernel::topods::Orientation;

/// OCCT TopBas_TestInterference — an interference with a parameter on a
/// curve, a boundary index and three orientations.
///
/// The OCCT default constructor is `= default` (members left uninitialized);
/// the Rust default fills them with zeros.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TestInterference {
    my_intersection: f64,
    my_boundary: i32,
    my_orientation: Orientation,
    my_transition: Orientation,
    my_b_transition: Orientation,
}

impl Default for TestInterference {
    fn default() -> Self {
        TestInterference {
            my_intersection: 0.0,
            my_boundary: 0,
            my_orientation: Orientation::Forward,
            my_transition: Orientation::Forward,
            my_b_transition: Orientation::Forward,
        }
    }
}

impl TestInterference {
    /// OCCT TopBas_TestInterference() — _0.cxx L30-31 (defaulted).
    pub fn new() -> Self {
        Self::default()
    }

    /// OCCT TopBas_TestInterference(Inters, Bound, Orient, Trans, BTrans) —
    /// _0.cxx L33-36.
    pub fn from_parts(
        inters: f64,
        bound: i32,
        orient: Orientation,
        trans: Orientation,
        b_trans: Orientation,
    ) -> Self {
        TestInterference {
            my_intersection: inters,
            my_boundary: bound,
            my_orientation: orient,
            my_transition: trans,
            my_b_transition: b_trans,
        }
    }

    /// OCCT Intersection(I) — hxx L76-79.
    pub fn set_intersection(&mut self, i: f64) {
        self.my_intersection = i;
    }

    /// OCCT Boundary(B) — hxx L81-84.
    pub fn set_boundary(&mut self, b: i32) {
        self.my_boundary = b;
    }

    /// OCCT Orientation(O) — hxx L86-89.
    pub fn set_orientation(&mut self, o: Orientation) {
        self.my_orientation = o;
    }

    /// OCCT Transition(Or) — hxx L91-94.
    pub fn set_transition(&mut self, or_: Orientation) {
        self.my_transition = or_;
    }

    /// OCCT BoundaryTransition(Or) — hxx L100-104.
    pub fn set_boundary_transition(&mut self, or_: Orientation) {
        self.my_b_transition = or_;
    }

    /// OCCT Intersection() — hxx L107-110.
    pub fn intersection(&self) -> f64 {
        self.my_intersection
    }

    /// OCCT ChangeIntersection() — hxx L112-115.
    pub fn change_intersection(&mut self) -> &mut f64 {
        &mut self.my_intersection
    }

    /// OCCT Boundary() — hxx L117-120.
    pub fn boundary(&self) -> i32 {
        self.my_boundary
    }

    /// OCCT ChangeBoundary() — hxx L122-125.
    pub fn change_boundary(&mut self) -> &mut i32 {
        &mut self.my_boundary
    }

    /// OCCT Orientation() — hxx L127-130.
    pub fn orientation(&self) -> Orientation {
        self.my_orientation
    }

    /// OCCT Transition() — hxx L132-135.
    pub fn transition(&self) -> Orientation {
        self.my_transition
    }

    /// OCCT BoundaryTransition() — hxx L137-140.
    pub fn boundary_transition(&self) -> Orientation {
        self.my_b_transition
    }
}
