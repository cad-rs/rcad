// OCCT TopClass_Classifier2d (TopClass_Classifier2d.pxx),
// BRepClass_FClass2dOfFClassifier, TopClass_FaceClassifier
// (TopClass_FaceClassifier.pxx) and BRepClass_FClassifier.
//
// The ray-casting 2D face classifier: cast a ray from the query point, count
// the crossings with the face's boundary edges, and derive the IN/OUT/ON
// state. Used as the precise fallback in FClass2d::Perform when the CSLib_2d
// polygon classification is uncertain (ON) or the wire is bad.

use glam::DVec2;
use rcad_kernel::geom::Line2d;

use crate::topalgo::shape_source::ShapeSource;
use crate::topalgo::brep_class::edge::ClassEdge;
use crate::topalgo::brep_class::face_explorer::FaceExplorer;
use crate::topalgo::brep_class::intersector::Intersector;
use crate::topalgo::brep_top_adaptor::fclass2d::State;
use crate::geomalgo::int_res2d::{IntersectionPoint, Position, TypeTrans};
use crate::geomalgo::top_trans::CurveTransition;

/// OCCT BRepClass_FClass2dOfFClassifier — the classifier used inside the
/// TopClass_FaceClassifier loop. State machine over the closest edge crossing.
pub struct FClass2dOfFClassifier {
    lin: Line2d,
    param: f64,
    tolerance: f64,
    state: State,
    first_compare: bool,
    first_trans: bool,
    closest: usize,
    is_set: bool,
    is_head_or_end: bool,
    trans: CurveTransition,
    intersector: Intersector,
}

impl FClass2dOfFClassifier {
    pub fn new() -> Self {
        FClass2dOfFClassifier {
            lin: Line2d::new(DVec2::ZERO, DVec2::X),
            param: 0.0,
            tolerance: 0.0,
            state: State::Unknown,
            first_compare: true,
            first_trans: true,
            closest: 0,
            is_set: false,
            is_head_or_end: false,
            trans: CurveTransition::new(),
            intersector: Intersector::new(),
        }
    }

    /// OCCT TopClass_Classifier2d::Reset (TopClass_Classifier2d.pxx L31-53).
    pub fn reset(&mut self, l: &Line2d, p: f64, tol: f64) {
        self.lin = l.clone();
        self.param = p;
        self.tolerance = tol;
        self.state = State::Unknown;
        self.first_compare = true;
        self.first_trans = true;
        self.closest = 0;
        self.is_set = true;
        self.is_head_or_end = false;
    }

    pub fn parameter(&self) -> f64 {
        self.param
    }

    pub fn state(&self) -> State {
        self.state
    }

    pub fn is_head_or_end(&self) -> bool {
        self.is_head_or_end
    }

    pub fn closest_intersection(&self) -> usize {
        self.closest
    }

    pub fn intersector(&self) -> &Intersector {
        &self.intersector
    }

    /// OCCT TopClass_Classifier2d::Compare (TopClass_Classifier2d.pxx L57-225).
    pub fn compare(&mut self, edge: &ClassEdge, or_: rcad_kernel::topods::Orientation, ds: &dyn ShapeSource) {
        // Intersect the edge and the segment.
        self.closest = 0;
        self.intersector.perform(&self.lin, self.param, self.tolerance, edge, ds);
        if !self.intersector.is_done() {
            return;
        }
        if self.intersector.nb_points() == 0 && self.intersector.nb_segments() == 0 {
            return;
        }

        // Find the closest point.
        let mut a_d_min = f64::MAX;
        let mut a_p_closest: Option<IntersectionPoint> = None;
        let a_nb_points = self.intersector.nb_points();
        for a_point in 1..=a_nb_points {
            let a_p_inter = self.intersector.point(a_point);
            // Test for ON.
            if a_p_inter.transition_of_first().position_on_curve() == Position::Head {
                self.closest = a_point;
                self.state = State::On;
                return;
            }
            let a_param_first = a_p_inter.param_on_first();
            if a_param_first < a_d_min {
                self.closest = a_point;
                a_p_closest = Some(a_p_inter.clone());
                a_d_min = a_param_first;
            }
        }

        // For the segments we only test the first point.
        let a_nb_segments = self.intersector.nb_segments();
        for a_segment in 1..=a_nb_segments {
            let a_seg_inter = self.intersector.segment(a_segment);
            let a_p_inter = a_seg_inter.first_point();
            if a_p_inter.transition_of_first().position_on_curve() == Position::Head {
                self.closest = a_nb_points + a_segment + a_segment - 1;
                self.state = State::On;
                return;
            }
            let a_param_first = a_p_inter.param_on_first();
            if a_param_first < a_d_min {
                self.closest = a_nb_points + a_segment + a_segment - 1;
                a_p_closest = Some(a_p_inter.clone());
                a_d_min = a_param_first;
            }
        }

        // If no point was found return.
        if self.closest == 0 {
            return;
        }

        // If the edge is INTERNAL or EXTERNAL.
        if or_ == rcad_kernel::topods::Orientation::Internal {
            self.state = State::In;
            return;
        } else if or_ == rcad_kernel::topods::Orientation::External {
            self.state = State::Out;
            return;
        }

        if !self.first_compare {
            if a_d_min > self.param {
                return;
            }
        }

        // Process the closest point aPClosest, found at aDMin on line.
        self.first_compare = false;
        if self.param > a_d_min {
            self.first_trans = true;
        }
        self.param = a_d_min;
        let a_p_closest = match a_p_closest {
            Some(p) => p,
            None => return,
        };
        let a_t2 = a_p_closest.transition_of_second();
        self.is_head_or_end = a_t2.position_on_curve() == Position::Head
            || a_t2.position_on_curve() == Position::End;

        // Transition on the segment.
        let mut a_seg_trans = rcad_kernel::topods::Orientation::Forward;
        let a_t1 = a_p_closest.transition_of_first();
        match a_t1.transition_type() {
            TypeTrans::In => {
                a_seg_trans = if or_ == rcad_kernel::topods::Orientation::Reversed {
                    rcad_kernel::topods::Orientation::Reversed
                } else {
                    rcad_kernel::topods::Orientation::Forward
                };
            }
            TypeTrans::Out => {
                a_seg_trans = if or_ == rcad_kernel::topods::Orientation::Reversed {
                    rcad_kernel::topods::Orientation::Forward
                } else {
                    rcad_kernel::topods::Orientation::Reversed
                };
            }
            TypeTrans::Touch => {
                match a_t1.situation() {
                    crate::geomalgo::int_res2d::Situation::Inside => {
                        a_seg_trans = if or_ == rcad_kernel::topods::Orientation::Reversed {
                            rcad_kernel::topods::Orientation::External
                        } else {
                            rcad_kernel::topods::Orientation::Internal
                        };
                    }
                    crate::geomalgo::int_res2d::Situation::Outside => {
                        a_seg_trans = if or_ == rcad_kernel::topods::Orientation::Reversed {
                            rcad_kernel::topods::Orientation::Internal
                        } else {
                            rcad_kernel::topods::Orientation::External
                        };
                    }
                    crate::geomalgo::int_res2d::Situation::Unknown => return,
                }
            }
            TypeTrans::Undecided => return,
        }

        if !self.is_head_or_end {
            // aPClosest is inside the edge.
            match a_seg_trans {
                rcad_kernel::topods::Orientation::Forward
                | rcad_kernel::topods::Orientation::External => self.state = State::Out,
                rcad_kernel::topods::Orientation::Reversed
                | rcad_kernel::topods::Orientation::Internal => self.state = State::In,
                _ => {}
            }
        } else {
            // aPClosest is Head or End of the edge: update the complex
            // transition.
            let (a_tang2d, a_norm2d, a_curv) =
                self.intersector.local_geometry(edge, ds, a_p_closest.param_on_second());
            if self.first_trans {
                self.trans.reset(self.lin.direction);
                self.first_trans = false;
            }
            let a_ort = if a_t2.position_on_curve() == Position::Head {
                rcad_kernel::topods::Orientation::Forward
            } else {
                rcad_kernel::topods::Orientation::Reversed
            };
            self.trans.compare(f64::EPSILON, a_tang2d, a_norm2d, a_curv, a_seg_trans, a_ort);
            self.state = self.trans.state_before();
        }
    }
}

/// OCCT TopClass_FaceClassifier + BRepClass_FClassifier — cast a ray from the
/// point and count crossings with the face boundary edges.
pub struct FClassifier {
    rejected: bool,
    nowires: bool,
    classifier: FClass2dOfFClassifier,
    // OCCT myEdge / myPosition / myEdgeParameter (kept for API completeness).
    edge: Option<ClassEdge>,
    position: Position,
    edge_parameter: f64,
}

impl FClassifier {
    pub fn new() -> Self {
        FClassifier {
            rejected: false,
            nowires: true,
            classifier: FClass2dOfFClassifier::new(),
            edge: None,
            position: Position::Middle,
            edge_parameter: 0.0,
        }
    }

    /// OCCT BRepClass_FClassifier::Perform(F, P, Tol) — classify a 2D point
    /// against the face using a BRepClass_FaceExplorer.
    pub fn perform(&mut self, ds: &dyn ShapeSource, face: usize, p: DVec2, tol: f64) {
        let mut fexp = FaceExplorer::new(ds, face);
        fexp.set_max_tolerance(0.1);
        fexp.set_use_bnd_box(false);

        // OCCT TopClass_FaceClassifier::Perform (TopClass_FaceClassifier.pxx
        // L37-144).
        let mut a_point = p;
        loop {
            a_point = fexp.check_point(ds, a_point);
            if a_point.is_finite() {
                break;
            }
        }
        self.rejected = fexp.reject(a_point);
        if self.rejected {
            return;
        }

        let mut a_line: Option<(Line2d, f64)> = fexp.segment(ds, a_point);
        self.nowires = true;

        while let Some((line, a_param)) = a_line {
            self.classifier.reset(&line, a_param, tol);

            fexp.init_wires();
            while fexp.more_wires() {
                self.nowires = false;
                if fexp.reject_wire(line.direction, self.classifier.parameter()) {
                    fexp.next_wire();
                    continue;
                }
                fexp.init_edges();
                while fexp.more_edges() {
                    if fexp.reject_edge(line.direction, self.classifier.parameter()) {
                        fexp.next_edge();
                        continue;
                    }
                    let an_edge_ori = fexp.current_edge_orientation();
                    if an_edge_ori == rcad_kernel::topods::Orientation::Forward
                        || an_edge_ori == rcad_kernel::topods::Orientation::Reversed
                    {
                        if let Some(an_edge) = fexp.current_edge(ds) {
                            self.classifier.compare(&an_edge, an_edge_ori, ds);
                            let a_closest_ind = self.classifier.closest_intersection();
                            if a_closest_ind != 0 {
                                let an_intersector = self.classifier.intersector();
                                let a_nb_pnts = an_intersector.nb_points();
                                self.edge = Some(an_edge);
                                if a_closest_ind <= a_nb_pnts {
                                    let a_p_inter = an_intersector.point(a_closest_ind);
                                    self.position = a_p_inter.transition_of_second().position_on_curve();
                                    self.edge_parameter = a_p_inter.param_on_second();
                                } else {
                                    let idx = a_closest_ind - a_nb_pnts;
                                    let a_p_inter = if idx & 1 == 1 {
                                        an_intersector
                                            .segment((idx + 1) / 2)
                                            .first_point()
                                            .clone()
                                    } else {
                                        an_intersector.segment((idx + 1) / 2).last_point().clone()
                                    };
                                    self.position = a_p_inter.transition_of_second().position_on_curve();
                                    self.edge_parameter = a_p_inter.param_on_second();
                                }
                            }
                        }
                    }
                    // If we are ON, we stop.
                    if self.classifier.state() == State::On {
                        return;
                    }
                    fexp.next_edge();
                }
                // If we are out of the wire we stop.
                if self.classifier.state() == State::Out {
                    return;
                }
                fexp.next_wire();
            }

            if !self.classifier.is_head_or_end() && self.classifier.state() != State::Unknown {
                break;
            }
            // Bad case for classification: try another segment.
            a_line = fexp.other_segment(ds, a_point);
        }
    }

    /// OCCT TopClass_FaceClassifier::State (TopClass_FaceClassifier.pxx L148-159).
    pub fn state(&self) -> State {
        if self.rejected {
            State::Out
        } else if self.nowires {
            State::In
        } else {
            self.classifier.state()
        }
    }

    /// OCCT BRepClass_FClassifier::Edge — the closest edge.
    pub fn edge(&self) -> Option<&ClassEdge> {
        self.edge.as_ref()
    }

    pub fn edge_parameter(&self) -> f64 {
        self.edge_parameter
    }

    pub fn position(&self) -> Position {
        self.position
    }
}
