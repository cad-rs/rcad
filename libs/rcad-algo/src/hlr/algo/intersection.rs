//! OCCT HLRAlgo_Intersection (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_Intersection.hxx` (L31-82) +
//! `HLRAlgo_Intersection.cxx` (L21-47) + `HLRAlgo_Intersection.lxx`
//! (L19-113).

use rcad_kernel::topods::{Orientation, State};

/// OCCT HLRAlgo_Intersection — an intersection on an edge to hide: a
/// parameter and a state (ON = on the face, OUT = above, IN = under).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intersection {
    my_orien: Orientation,
    my_seg_index: i32,
    my_index: i32,
    my_level: i32,
    my_param: f64,
    my_toler: f32,
    my_state: State,
}

impl Default for Intersection {
    fn default() -> Self {
        // OCCT default ctor (cxx L21-28): myOrien left uninitialized
        // (TopAbs default 0 = Forward in the Rust enum), myState likewise.
        Intersection {
            my_orien: Orientation::Forward,
            my_seg_index: 0,
            my_index: 0,
            my_level: 0,
            my_param: 0.0,
            my_toler: 0.0,
            my_state: State::In,
        }
    }
}

impl Intersection {
    /// OCCT HLRAlgo_Intersection() — cxx L21-28.
    pub fn new() -> Self {
        Self::default()
    }

    /// OCCT HLRAlgo_Intersection(Ori, Lev, SegInd, Ind, P, Tol, S) —
    /// cxx L32-47.
    #[allow(clippy::too_many_arguments)]
    pub fn from_parts(
        ori: Orientation,
        lev: i32,
        seg_ind: i32,
        ind: i32,
        p: f64,
        tol: f32,
        s: State,
    ) -> Self {
        Intersection {
            my_orien: ori,
            my_seg_index: seg_ind,
            my_index: ind,
            my_level: lev,
            my_param: p,
            my_toler: tol,
            my_state: s,
        }
    }

    /// OCCT Orientation(Ori) — lxx L19-22.
    pub fn set_orientation(&mut self, ori: Orientation) {
        self.my_orien = ori;
    }

    /// OCCT Orientation() — lxx L26-29.
    pub fn orientation(&self) -> Orientation {
        self.my_orien
    }

    /// OCCT Level(Lev) — lxx L33-36.
    pub fn set_level(&mut self, lev: i32) {
        self.my_level = lev;
    }

    /// OCCT Level() — lxx L40-43.
    pub fn level(&self) -> i32 {
        self.my_level
    }

    /// OCCT SegIndex(SegInd) — lxx L47-50.
    pub fn set_seg_index(&mut self, seg_ind: i32) {
        self.my_seg_index = seg_ind;
    }

    /// OCCT SegIndex() — lxx L54-57.
    pub fn seg_index(&self) -> i32 {
        self.my_seg_index
    }

    /// OCCT Index(Ind) — lxx L61-64.
    pub fn set_index(&mut self, ind: i32) {
        self.my_index = ind;
    }

    /// OCCT Index() — lxx L68-71.
    pub fn index(&self) -> i32 {
        self.my_index
    }

    /// OCCT Parameter(P) — lxx L75-78.
    pub fn set_parameter(&mut self, p: f64) {
        self.my_param = p;
    }

    /// OCCT Parameter() — lxx L82-85.
    pub fn parameter(&self) -> f64 {
        self.my_param
    }

    /// OCCT Tolerance(T) — lxx L89-92.
    pub fn set_tolerance(&mut self, t: f32) {
        self.my_toler = t;
    }

    /// OCCT Tolerance() — lxx L96-99.
    pub fn tolerance(&self) -> f32 {
        self.my_toler
    }

    /// OCCT State(S) — lxx L103-106.
    pub fn set_state(&mut self, s: State) {
        self.my_state = s;
    }

    /// OCCT State() — lxx L110-113.
    pub fn state(&self) -> State {
        self.my_state
    }
}
