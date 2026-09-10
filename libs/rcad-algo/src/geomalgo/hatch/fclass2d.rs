//! OCCT Geom2dHatch_FClass2dOfClassifier (TKGeomAlgo/Geom2dHatch).
//!
//! Geom2dHatch_FClass2dOfClassifier.hxx L34-85 + .cxx L27-74 — the classifier
//! fed by the TopClass_FaceClassifier loop for the 2D hatching.  Its Reset /
//! Compare are the TopClass_Classifier2d.pxx template bodies (L31-53 /
//! L57-225) instantiated with TheEdge = Geom2dAdaptor_Curve (rcad `Curve2d`)
//! and TheIntersector = Geom2dHatch_Intersector.  (Translation form follows
//! the landed BRepClass precedent: the pxx bodies are inlined into the class
//! methods.)

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Line2d};
use rcad_kernel::topods::{Orientation, State};

use crate::geomalgo::hatch::intersector::HatchIntersector;
use crate::geomalgo::int_res2d::{
    IntersectionPoint as IntRes2dIntersectionPoint, Position, Situation, TypeTrans,
};
use crate::geomalgo::top_trans::CurveTransition;

/// OCCT Geom2dHatch_FClass2dOfClassifier — the per-edge classification state
/// machine over the closest crossing of the probing segment.
pub struct FClass2dOfClassifier {
    /// OCCT bool myIsSet.
    my_is_set: bool,
    /// OCCT bool myFirstCompare.
    my_first_compare: bool,
    /// OCCT bool myFirstTrans.
    my_first_trans: bool,
    /// OCCT gp_Lin2d myLin.
    my_lin: Line2d,
    /// OCCT double myParam.
    my_param: f64,
    /// OCCT double myTolerance.
    my_tolerance: f64,
    /// OCCT TopTrans_CurveTransition myTrans.
    my_trans: CurveTransition<DVec2>,
    /// OCCT Geom2dHatch_Intersector myIntersector.
    my_intersector: HatchIntersector,
    /// OCCT int myClosest.
    my_closest: usize,
    /// OCCT TopAbs_State myState.
    my_state: State,
    /// OCCT bool myIsHeadOrEnd.
    my_is_head_or_end: bool,
}

impl FClass2dOfClassifier {
    /// OCCT Geom2dHatch_FClass2dOfClassifier() (cxx L27-37) — creates an
    /// undefined classifier.
    pub fn new() -> Self {
        FClass2dOfClassifier {
            my_is_set: false,
            my_first_compare: true,
            my_first_trans: true,
            // OCCT gp_Lin2d default-constructs on the X axis.
            my_lin: Line2d::new(DVec2::ZERO, DVec2::X),
            my_param: 0.0,
            my_tolerance: 0.0,
            my_trans: CurveTransition::new(),
            my_intersector: HatchIntersector::new(),
            my_closest: 0,
            my_state: State::Unknown,
            my_is_head_or_end: false,
        }
    }

    /// OCCT Reset (cxx L41-55) — the TopClass_Classifier2d::Reset body
    /// (TopClass_Classifier2d.pxx L31-53).  Starts a classification process:
    /// the point to classify is the origin of the line <L>; <P> is the
    /// original length of the segment on <L>; <Tol> is the tolerance.
    pub fn reset(&mut self, l: &Line2d, p: f64, tol: f64) {
        self.my_lin = *l;
        self.my_param = p;
        self.my_tolerance = tol;
        self.my_state = State::Unknown;
        self.my_first_compare = true;
        self.my_first_trans = true;
        self.my_closest = 0;
        self.my_is_set = true;
        self.my_is_head_or_end = false;
    }

    /// OCCT Compare (cxx L59-74) — the TopClass_Classifier2d::Compare body
    /// (TopClass_Classifier2d.pxx L57-225).  Updates the classification
    /// process with the edge <E> from the boundary.
    // OCCT initializes aSegTrans to TopAbs_FORWARD before the switch; some
    // switch arms return before it is read (form-preserved initial value).
    #[allow(unused_assignments)]
    pub fn compare(&mut self, the_edge: &Curve2d, the_or: Orientation) {
        // Intersect the edge and the segment.
        self.my_closest = 0;
        self.my_intersector
            .perform(&self.my_lin, self.my_param, self.my_tolerance, the_edge);
        if !self.my_intersector.is_done() {
            return;
        }
        if (self.my_intersector.nb_points() == 0) && (self.my_intersector.nb_segments() == 0) {
            return;
        }

        // Find the closest point.
        let mut a_p_closest: Option<IntRes2dIntersectionPoint> = None;

        let mut a_d_min = f64::MAX; // RealLast()
        let a_nb_points = self.my_intersector.nb_points();
        for a_point in 1..=a_nb_points {
            let a_p_inter = self.my_intersector.point(a_point);
            // Test for ON.
            if a_p_inter.transition_of_first().position_on_curve() == Position::Head {
                self.my_closest = a_point;
                self.my_state = State::On;
                return;
            }
            let a_param_first = a_p_inter.param_on_first();
            if a_param_first < a_d_min {
                self.my_closest = a_point;
                a_p_closest = Some(a_p_inter.clone());
                a_d_min = a_param_first;
            }
        }

        // For the segments we only test the first point.
        let a_nb_segments = self.my_intersector.nb_segments();
        for a_segment in 1..=a_nb_segments {
            let a_seg_inter = self.my_intersector.segment(a_segment);
            let a_p_inter = a_seg_inter.first_point();
            if a_p_inter.transition_of_first().position_on_curve() == Position::Head {
                self.my_closest = a_nb_points + a_segment + a_segment - 1;
                self.my_state = State::On;
                return;
            }
            let a_param_first = a_p_inter.param_on_first();
            if a_param_first < a_d_min {
                self.my_closest = a_nb_points + a_segment + a_segment - 1;
                a_p_closest = Some(a_p_inter.clone());
                a_d_min = a_param_first;
            }
        }

        // If no point was found return.
        if self.my_closest == 0 {
            return;
        }

        // If the edge is INTERNAL or EXTERNAL, no problem.
        if the_or == Orientation::Internal {
            self.my_state = State::In;
            return;
        } else if the_or == Orientation::External {
            self.my_state = State::Out;
            return;
        }

        if !self.my_first_compare {
            if a_d_min > self.my_param {
                return;
            }
        }

        // Process the closest point aPClosest, found at aDMin on line.
        self.my_first_compare = false;

        if self.my_param > a_d_min {
            self.my_first_trans = true;
        }

        self.my_param = a_d_min;
        let a_p_closest = match a_p_closest {
            Some(p) => p,
            // OCCT dereferences the pointer unconditionally; a closest point
            // is always registered when myClosest != 0.
            None => return,
        };
        let a_t2_position = a_p_closest.transition_of_second().position_on_curve();
        self.my_is_head_or_end =
            (a_t2_position == Position::Head) || (a_t2_position == Position::End);

        // Transition on the segment.
        let mut a_seg_trans = Orientation::Forward;

        let a_t1 = a_p_closest.transition_of_first();
        match a_t1.transition_type() {
            TypeTrans::In => {
                a_seg_trans = if the_or == Orientation::Reversed {
                    Orientation::Reversed
                } else {
                    Orientation::Forward
                };
            }
            TypeTrans::Out => {
                a_seg_trans = if the_or == Orientation::Reversed {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
            }
            TypeTrans::Touch => match a_t1.situation() {
                Situation::Inside => {
                    a_seg_trans = if the_or == Orientation::Reversed {
                        Orientation::External
                    } else {
                        Orientation::Internal
                    };
                }
                Situation::Outside => {
                    a_seg_trans = if the_or == Orientation::Reversed {
                        Orientation::Internal
                    } else {
                        Orientation::External
                    };
                }
                Situation::Unknown => return,
            },
            TypeTrans::Undecided => return,
        }

        // Are we inside the edge?
        if !self.my_is_head_or_end {
            // aPClosest is inside the edge.
            match a_seg_trans {
                Orientation::Forward | Orientation::External => self.my_state = State::Out,
                Orientation::Reversed | Orientation::Internal => self.my_state = State::In,
            }
        } else {
            // aPClosest is Head or End of the edge: update the complex
            // transition.
            // OCCT theIntersector.LocalGeometry(theEdge,
            // aPClosest->ParamOnSecond(), aTang2d, aNorm2d, aCurv).
            let (a_tang2d, a_norm2d, a_curv) = self
                .my_intersector
                .local_geometry(the_edge, a_p_closest.param_on_second());
            if self.my_first_trans {
                // OCCT builds 3D gp_Dir from the line direction (Z = 0); the
                // rcad TopTrans_CurveTransition is carried in 2D.
                self.my_trans.reset(self.my_lin.direction);
                self.my_first_trans = false;
            }

            let a_ort = if a_t2_position == Position::Head {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            // OCCT theTrans.Compare(RealEpsilon(), aTang, aNorm, aCurv,
            //                       aSegTrans, aOrt).
            self.my_trans
                .compare(f64::EPSILON, a_tang2d, a_norm2d, a_curv, a_seg_trans, a_ort);
            self.my_state = self.my_trans.state_before();
        }
    }

    /// OCCT Parameter (hxx L54) — the current value of the parameter.
    pub fn parameter(&self) -> f64 {
        self.my_param
    }

    /// OCCT Intersector (hxx L57) — the intersecting algorithm.
    pub fn intersector(&self) -> &HatchIntersector {
        &self.my_intersector
    }

    /// OCCT Intersector& (hxx L57) — mutable access.
    pub fn intersector_mut(&mut self) -> &mut HatchIntersector {
        &mut self.my_intersector
    }

    /// OCCT ClosestIntersection (hxx L63) — 0 if the last compared edge had
    /// no relevant intersection, else the index of this intersection.
    pub fn closest_intersection(&self) -> usize {
        self.my_closest
    }

    /// OCCT State (hxx L66) — the current state of the point.
    pub fn state(&self) -> State {
        self.my_state
    }

    /// OCCT IsHeadOrEnd (hxx L71) — true when the closest intersection point
    /// represents the head or end of the edge.
    pub fn is_head_or_end(&self) -> bool {
        self.my_is_head_or_end
    }
}

impl Default for FClass2dOfClassifier {
    fn default() -> Self {
        Self::new()
    }
}
