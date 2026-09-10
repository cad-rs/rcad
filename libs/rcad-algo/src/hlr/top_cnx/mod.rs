//! OCCT TopCnx package (TKHLR) — cumulated transition on an edge.
//!
//! 1:1 translation of `TopCnx_EdgeFaceTransition.hxx` (L17-71) +
//! `TopCnx_EdgeFaceTransition.cxx` (L17-134).  `TopTrans_CurveTransition`
//! maps to the shared `geomalgo::top_trans::CurveTransition<DVec3>` port.

use glam::DVec3;
use rcad_kernel::topods::{Orientation, State};

use crate::geomalgo::top_trans::CurveTransition;

/// OCCT TopCnx_EdgeFaceTransition — computes the cumulated transition for
/// interferences on an edge.
pub struct EdgeFaceTransition {
    my_curve_transition: CurveTransition<DVec3>,
    nb_bound_forward: i32,
    nb_bound_reversed: i32,
}

impl Default for EdgeFaceTransition {
    fn default() -> Self {
        Self::new()
    }
}

impl EdgeFaceTransition {
    /// OCCT TopCnx_EdgeFaceTransition() — TopCnx_EdgeFaceTransition.cxx
    /// L20-25.
    pub fn new() -> Self {
        EdgeFaceTransition {
            my_curve_transition: CurveTransition::new(),
            nb_bound_forward: 0,
            nb_bound_reversed: 0,
        }
    }

    /// OCCT Reset(Tgt, Norm, Curv) — TopCnx_EdgeFaceTransition.cxx L28-33.
    pub fn reset(&mut self, tgt: DVec3, norm: DVec3, curv: f64) {
        self.my_curve_transition.reset_3d(tgt, norm, curv);
        self.nb_bound_forward = 0;
        self.nb_bound_reversed = 0;
    }

    /// OCCT Reset(Tgt) — TopCnx_EdgeFaceTransition.cxx L36-41 (a linear
    /// edge).
    pub fn reset_linear(&mut self, tgt: DVec3) {
        self.my_curve_transition.reset(tgt);
        self.nb_bound_forward = 0;
        self.nb_bound_reversed = 0;
    }

    /// OCCT AddInterference(Tole, Tang, Norm, Curv, Or, Tr, BTr) —
    /// TopCnx_EdgeFaceTransition.cxx L44-70.
    pub fn add_interference(
        &mut self,
        tole: f64,
        tang: DVec3,
        norm: DVec3,
        curv: f64,
        or_: Orientation,
        tr: Orientation,
        b_tr: Orientation,
    ) {
        self.my_curve_transition
            .compare(tole, tang, norm, curv, tr, or_);
        match b_tr {
            Orientation::Forward => {
                self.nb_bound_forward += 1;
            }
            Orientation::Reversed => {
                self.nb_bound_reversed += 1;
            }
            Orientation::Internal | Orientation::External => {}
        }
    }

    /// OCCT Transition() — TopCnx_EdgeFaceTransition.cxx L73-124: the
    /// cumulated transition from the states before / after the interferences.
    pub fn transition(&self) -> Orientation {
        let bef = self.my_curve_transition.state_before();
        let aft = self.my_curve_transition.state_after();
        if bef == State::In {
            if aft == State::In {
                return Orientation::Internal;
            } else if aft == State::Out {
                return Orientation::Reversed;
            }
        } else if bef == State::Out {
            if aft == State::In {
                return Orientation::Forward;
            } else if aft == State::Out {
                return Orientation::External;
            }
        }
        Orientation::Internal
    }

    /// OCCT BoundaryTransition() — TopCnx_EdgeFaceTransition.cxx L127-134.
    /// The last two branches return the same value in OCCT (kept verbatim).
    pub fn boundary_transition(&self) -> Orientation {
        if self.nb_bound_forward > self.nb_bound_reversed {
            Orientation::Forward
        } else if self.nb_bound_forward < self.nb_bound_reversed {
            Orientation::Reversed
        } else if (self.nb_bound_reversed % 2) == 0 {
            Orientation::External
        } else {
            Orientation::External
        }
    }
}
