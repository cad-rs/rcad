// OCCT TopTrans_CurveTransition (TopTrans_CurveTransition.cxx / .hxx)
// Complex transition of a curve relative to a boundary near an interference.
//
// Used by TopClass_Classifier2d::Compare (the FClass2d face-classifier
// fallback) to combine the transitions at the head/end of a boundary edge
// into a single IN/OUT state.

pub mod surface_transition; // TopTrans_SurfaceTransition (surface variant)

use glam::{DVec2, DVec3};
use rcad_kernel::topods::{Orientation, State};

const GREATER: i32 = 1;
const SAME: i32 = 0;

/// Minimal vector interface for the transition vectors (OCCT gp_Dir is 3D;
/// the FClass2d face classifier uses the 2D specialization).
pub trait TransitionVec: Copy {
    const ZERO: Self;
    fn dot_v(self, other: Self) -> f64;
    fn neg_v(self) -> Self;
}

impl TransitionVec for DVec2 {
    const ZERO: Self = DVec2::ZERO;
    fn dot_v(self, other: Self) -> f64 {
        self.dot(other)
    }
    fn neg_v(self) -> Self {
        -self
    }
}

impl TransitionVec for DVec3 {
    const ZERO: Self = DVec3::ZERO;
    fn dot_v(self, other: Self) -> f64 {
        self.dot(other)
    }
    fn neg_v(self) -> Self {
        -self
    }
}
const LOWER: i32 = -1;

/// OCCT TopAbs::Reverse — reverses an orientation.
fn reverse_orientation(or: Orientation) -> Orientation {
    match or {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}

/// OCCT TopTrans_CurveTransition — accumulates the transitions of several
/// interference points on an intersected curve and gives the state before /
/// after the complex interference.
pub struct CurveTransition<T: TransitionVec> {
    // Reference curve elements (the intersecting straight line).
    my_tgt: T,
    my_norm: T,
    my_curv: f64,
    // Init = first comparison flag.
    init: bool,
    // First (last) interference elements.
    tgt_first: T,
    norm_first: T,
    curv_first: f64,
    tran_first: Orientation,
    tgt_last: T,
    norm_last: T,
    curv_last: f64,
    tran_last: Orientation,
}

impl<T: TransitionVec> CurveTransition<T> {
    pub fn new() -> Self {
        CurveTransition {
            my_tgt: T::ZERO,
            my_norm: T::ZERO,
            my_curv: 0.0,
            init: false,
            tgt_first: T::ZERO,
            norm_first: T::ZERO,
            curv_first: 0.0,
            tran_first: Orientation::Forward,
            tgt_last: T::ZERO,
            norm_last: T::ZERO,
            curv_last: 0.0,
            tran_last: Orientation::Forward,
        }
    }

    /// OCCT Reset(Tgt) — initializer with the intersecting straight line.
    pub fn reset(&mut self, tgt: T) {
        self.my_tgt = tgt;
        self.my_curv = 0.0;
        self.init = true;
    }

    /// OCCT Reset(Tgt, Norm, Curv) (L42-48) — initializer with the elements of
    /// the intersecting curve.
    pub fn reset_3d(&mut self, tgt: T, norm: T, curv: f64) {
        self.my_tgt = tgt;
        self.my_norm = norm;
        self.my_curv = curv;
        self.init = true;
    }

    /// OCCT Compare (TopTrans_CurveTransition.cxx L69-273).
    /// `st` is the segment transition (how the curve crosses the boundary),
    /// `or_` the orientation of the interference on the boundary edge.
    pub fn compare(
        &mut self,
        tole: f64,
        t: T,
        n: T,
        c: f64,
        st: Orientation,
        or_: Orientation,
    ) {
        // S is the transition, O the orientation of the intersection on the
        // boundary.
        let mut s = st;
        let mut o = or_;

        // Adjustment for INTERNAL transition.
        if s == Orientation::Internal {
            if t.dot_v(self.my_tgt) < 0.0 {
                s = reverse_orientation(o);
            } else {
                s = o;
            }
        }

        if self.init {
            // First comparison for this complex transition.
            self.init = false;
            self.tgt_first = t;
            self.norm_first = n;
            self.curv_first = c;
            self.tran_first = s;
            self.tgt_last = t;
            self.norm_last = n;
            self.curv_last = c;
            self.tran_last = s;
            match o {
                // Interference at the end of the edge: reverse the tangent.
                Orientation::Reversed => {
                    self.tgt_first = self.tgt_first.neg_v();
                    self.tgt_last = self.tgt_last.neg_v();
                }
                // Interference in the middle of the edge: reverse depending on
                // the position of the reference tangent.
                Orientation::Internal => {
                    if self.my_tgt.dot_v(t) > 0.0 {
                        self.tgt_first = self.tgt_first.neg_v();
                    } else {
                        self.tgt_last = self.tgt_last.neg_v();
                    }
                }
                Orientation::Forward | Orientation::External => {}
            }
        } else {
            // Compare with the existing first and last transition.
            let mut first_set = false;
            let mut cos_ang_with_t = self.my_tgt.dot_v(t);
            match o {
                Orientation::Reversed => cos_ang_with_t = -cos_ang_with_t,
                Orientation::Internal => {
                    if cos_ang_with_t > 0.0 {
                        cos_ang_with_t = -cos_ang_with_t;
                    }
                }
                Orientation::Forward | Orientation::External => {}
            }
            let cos_ang_with_1 = self.my_tgt.dot_v(self.tgt_first);

            match Self::compare_angles(cos_ang_with_t, cos_ang_with_1, tole) {
                LOWER => {
                    // The angle is greater than the first: the new one becomes
                    // the first.
                    first_set = true;
                    self.tgt_first = t;
                    match o {
                        Orientation::Reversed => {
                            self.tgt_first = self.tgt_first.neg_v()
                        }
                        Orientation::Internal => {
                            if self.my_tgt.dot_v(t) > 0.0 {
                                self.tgt_first = self.tgt_first.neg_v();
                            }
                        }
                        Orientation::Forward | Orientation::External => {}
                    }
                    self.norm_first = n;
                    self.curv_first = c;
                    self.tran_first = s;
                }
                SAME => {
                    // Same angle: look at the curvature.
                    if self.is_before(tole, cos_ang_with_t, n, c, self.norm_first, self.curv_first)
                    {
                        first_set = true;
                        self.tgt_first = t;
                        match o {
                            Orientation::Reversed => {
                            self.tgt_first = self.tgt_first.neg_v()
                        }
                            Orientation::Internal => {
                                if self.my_tgt.dot_v(t) > 0.0 {
                                    self.tgt_first = self.tgt_first.neg_v();
                                }
                            }
                            Orientation::Forward | Orientation::External => {}
                        }
                        self.norm_first = n;
                        self.curv_first = c;
                        self.tran_first = s;
                    }
                }
                GREATER => {}
                _ => unreachable!(),
            }

            if !first_set || o == Orientation::Internal {
                // In tangency cases the first can also be the last.
                if o == Orientation::Internal {
                    cos_ang_with_t = -cos_ang_with_t;
                }
                let cos_ang_with_2 = self.my_tgt.dot_v(self.tgt_last);

                match Self::compare_angles(cos_ang_with_t, cos_ang_with_2, tole) {
                    GREATER => {
                        // The angle is lower than the last: the new one becomes
                        // the last.
                        self.tgt_last = t;
                        match o {
                            Orientation::Reversed => {
                            self.tgt_last = self.tgt_last.neg_v()
                        }
                            Orientation::Internal => {
                                if self.my_tgt.dot_v(t) < 0.0 {
                                    self.tgt_last = self.tgt_last.neg_v();
                                }
                            }
                            Orientation::Forward | Orientation::External => {}
                        }
                        self.norm_last = n;
                        self.curv_last = c;
                        self.tran_last = s;
                    }
                    SAME => {
                        // Same angle: look at the curvature.
                        if self.is_before(
                            tole,
                            cos_ang_with_t,
                            self.norm_last,
                            self.curv_last,
                            n,
                            c,
                        ) {
                            self.tgt_last = t;
                            match o {
                                Orientation::Reversed => {
                            self.tgt_last = self.tgt_last.neg_v()
                        }
                                Orientation::Internal => {
                                    if self.my_tgt.dot_v(t) < 0.0 {
                                        self.tgt_last = self.tgt_last.neg_v();
                                    }
                                }
                                Orientation::Forward | Orientation::External => {}
                            }
                            self.norm_last = n;
                            self.curv_last = c;
                            self.tran_last = s;
                        }
                    }
                    LOWER => {}
                    _ => unreachable!(),
                }
            }
        }
    }

    /// OCCT StateBefore — state of the curve before the interference.
    pub fn state_before(&self) -> State {
        if self.init {
            return State::Unknown;
        }
        match self.tran_first {
            Orientation::Forward | Orientation::External => State::Out,
            Orientation::Reversed | Orientation::Internal => State::In,
            _ => State::Out,
        }
    }

    /// OCCT StateAfter — state of the curve after the interference.
    pub fn state_after(&self) -> State {
        if self.init {
            return State::Unknown;
        }
        match self.tran_last {
            Orientation::Forward | Orientation::Internal => State::In,
            Orientation::Reversed | Orientation::External => State::Out,
            _ => State::Out,
        }
    }

    /// OCCT IsBefore (TopTrans_CurveTransition.cxx L327-424) — true if the
    /// interference (T1,C1) happens before (T2,C2) in the crossing order.
    fn is_before(
        &self,
        tole: f64,
        cos_angl: f64,
        n1: T,
        c1: f64,
        n2: T,
        c2: f64,
    ) -> bool {
        let tn1 = self.my_tgt.dot_v(n1);
        let tn2 = self.my_tgt.dot_v(n2);
        let mut one_before = false;

        if tn1.abs() <= tole || tn2.abs() <= tole {
            // Tangent: the first is the interference with the curvature
            // nearest to the reference.
            if self.my_curv == 0.0 {
                // The reference is straight; the first has the lowest curvature.
                if c1 < c2 {
                    one_before = true;
                }
                if cos_angl > 0.0 {
                    one_before = !one_before;
                }
            } else {
                // The reference is curved; the first has the nearest curvature
                // in the direction.
                let delta_c1 = if c1 == 0.0 || self.my_curv == 0.0 {
                    c1 - self.my_curv
                } else {
                    (c1 - self.my_curv) * (n1.dot_v(self.my_norm))
                };
                let delta_c2 = if c2 == 0.0 || self.my_curv == 0.0 {
                    c2 - self.my_curv
                } else {
                    (c2 - self.my_curv) * (n2.dot_v(self.my_norm))
                };
                if delta_c1 < delta_c2 {
                    one_before = true;
                }
                if cos_angl > 0.0 {
                    one_before = !one_before;
                }
            }
        } else if tn1 < 0.0 {
            // Before the first interference we are inside the curvature.
            if tn2 > 0.0 {
                // Before the second we are outside the curvature.
                // The first interference is before.  /* ->)( */
                one_before = true;
            } else if c1 > c2 {
                // Both inside; choose the greater curvature. /* ->)) */
                one_before = true;
            }
        } else if tn1 > 0.0 {
            // Before the first interference we are outside the curvature.
            if tn2 > 0.0 {
                // Before the second we are outside the curvature. /* ->(( */
                if c1 < c2 {
                    one_before = true;
                }
            }
        }
        one_before
    }

    /// OCCT Compare(Ang1, Ang2, Tole) — compare two cosines with tolerance.
    fn compare_angles(ang1: f64, ang2: f64, tole: f64) -> i32 {
        let mut res = SAME;
        if ang1 - ang2 > tole {
            res = GREATER;
        } else if ang2 - ang1 > tole {
            res = LOWER;
        }
        res
    }
}
