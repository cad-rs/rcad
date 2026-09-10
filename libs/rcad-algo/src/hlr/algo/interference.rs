//! OCCT HLRAlgo_Interference (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_Interference.hxx` (L29-158, fully inline) +
//! `HLRAlgo_Interference_0.cxx` (L24-39).

use rcad_kernel::topods::Orientation;

use super::coincidence::Coincidence;
use super::intersection::Intersection;

/// OCCT HLRAlgo_Interference — an intersection on the hidden edge plus the
/// boundary (coincidence) data and the three orientations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interference {
    my_intersection: Intersection,
    my_boundary: Coincidence,
    my_orientation: Orientation,
    my_transition: Orientation,
    my_b_transition: Orientation,
}

impl Default for Interference {
    fn default() -> Self {
        Interference {
            my_intersection: Intersection::new(),
            my_boundary: Coincidence::new(),
            my_orientation: Orientation::Forward,
            my_transition: Orientation::Forward,
            my_b_transition: Orientation::Forward,
        }
    }
}

impl Interference {
    /// OCCT HLRAlgo_Interference() — _0.cxx L24.
    pub fn new() -> Self {
        Self::default()
    }

    /// OCCT HLRAlgo_Interference(Inters, Bound, Orient, Trans, BTrans) —
    /// _0.cxx L28-39.
    pub fn from_parts(
        inters: Intersection,
        bound: Coincidence,
        orient: Orientation,
        trans: Orientation,
        b_trans: Orientation,
    ) -> Self {
        Interference {
            my_intersection: inters,
            my_boundary: bound,
            my_orientation: orient,
            my_transition: trans,
            my_b_transition: b_trans,
        }
    }

    /// OCCT Intersection(I) — hxx L78-81.
    pub fn set_intersection(&mut self, i: Intersection) {
        self.my_intersection = i;
    }

    /// OCCT Boundary(B) — hxx L85-88.
    pub fn set_boundary(&mut self, b: Coincidence) {
        self.my_boundary = b;
    }

    /// OCCT Orientation(O) — hxx L92-95.
    pub fn set_orientation(&mut self, o: Orientation) {
        self.my_orientation = o;
    }

    /// OCCT Transition(Or) — hxx L99-102.
    pub fn set_transition(&mut self, or_: Orientation) {
        self.my_transition = or_;
    }

    /// OCCT BoundaryTransition(Or) — hxx L106-109.
    pub fn set_boundary_transition(&mut self, or_: Orientation) {
        self.my_b_transition = or_;
    }

    /// OCCT Intersection() — hxx L113-116.
    pub fn intersection(&self) -> &Intersection {
        &self.my_intersection
    }

    /// OCCT ChangeIntersection() — hxx L120-123.
    pub fn change_intersection(&mut self) -> &mut Intersection {
        &mut self.my_intersection
    }

    /// OCCT Boundary() — hxx L127-130.
    pub fn boundary(&self) -> &Coincidence {
        &self.my_boundary
    }

    /// OCCT ChangeBoundary() — hxx L134-137.
    pub fn change_boundary(&mut self) -> &mut Coincidence {
        &mut self.my_boundary
    }

    /// OCCT Orientation() — hxx L141-144.
    pub fn orientation(&self) -> Orientation {
        self.my_orientation
    }

    /// OCCT Transition() — hxx L148-151.
    pub fn transition(&self) -> Orientation {
        self.my_transition
    }

    /// OCCT BoundaryTransition() — hxx L155-158.
    pub fn boundary_transition(&self) -> Orientation {
        self.my_b_transition
    }
}
