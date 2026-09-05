//! OCCT Geom2dHatch_Classifier (TKGeomAlgo/Geom2dHatch).
//!
//! Geom2dHatch_Classifier.hxx L34-83 + .cxx L29-82 — classifies a 2D point
//! against the face described by Geom2dHatch_Elements.  Its Perform / State /
//! Edge / EdgeParameter are the TopClass_FaceClassifier.pxx template bodies
//! (L29-144 / L148-159 / L163-168 / L172-176) instantiated with
//! TheFaceExplorer = Geom2dHatch_Elements, TheEdge = Geom2dAdaptor_Curve
//! (rcad `Curve2d`) and FClass2d = Geom2dHatch_FClass2dOfClassifier.
//! (Translation form follows the landed BRepClass precedent: the pxx bodies
//! are inlined into the class methods.)

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Line2d};
use rcad_kernel::topods::{Orientation, State};

use crate::geomalgo::hatch::elements::HatchElements;
use crate::geomalgo::hatch::fclass2d::FClass2dOfClassifier;
use crate::geomalgo::int_res2d::Position;

/// OCCT Geom2dHatch_Classifier — the 2D point classifier on a hatched face.
pub struct Classifier {
    /// OCCT Geom2dHatch_FClass2dOfClassifier myClassifier.
    my_classifier: FClass2dOfClassifier,
    /// OCCT Geom2dAdaptor_Curve myEdge — default-constructed until Perform
    /// stores the closest edge (rcad: Option).
    my_edge: Option<Curve2d>,
    /// OCCT double myEdgeParameter.
    my_edge_parameter: f64,
    /// OCCT IntRes2d_Position myPosition — left uninitialized by the OCCT
    /// constructor; rcad starts at a neutral Middle.
    my_position: Position,
    /// OCCT bool rejected.
    rejected: bool,
    /// OCCT bool nowires.
    nowires: bool,
}

impl Classifier {
    /// OCCT Geom2dHatch_Classifier() (cxx L29-34) — empty constructor,
    /// undefined algorithm.
    pub fn new() -> Self {
        Classifier {
            my_classifier: FClass2dOfClassifier::new(),
            my_edge: None,
            my_edge_parameter: 0.0,
            my_position: Position::Middle,
            rejected: false,
            nowires: true,
        }
    }

    /// OCCT Geom2dHatch_Classifier(F, P, Tol) (cxx L38-46) — creates an
    /// algorithm to classify the point P with tolerance <Tol> on the face
    /// described by <F>.
    pub fn with_perform(f: &mut HatchElements, p: DVec2, tol: f64) -> Self {
        let mut r = Classifier::new();
        r.perform(f, p, tol);
        r
    }

    /// OCCT Perform (cxx L50-61) — the TopClass_FaceClassifier::Perform body
    /// (TopClass_FaceClassifier.pxx L29-144).  Classifies the point P with
    /// tolerance <Tol> on the face described by <F>.
    // OCCT declares anEdge / anEdgeOri / aPInter once and assigns them inside
    // the loops; the initial values stay form-preserved (hence the allow).
    #[allow(unused_assignments)]
    pub fn perform(&mut self, the_fexp: &mut HatchElements, the_p: DVec2, the_tol: f64) {
        let mut a_point = the_p;
        let mut a_res_of_point_check = false;
        while a_res_of_point_check == false {
            a_res_of_point_check = the_fexp.check_point(&mut a_point);
        }

        // Test for rejection.
        self.rejected = the_fexp.reject(a_point);

        if self.rejected {
            return;
        }

        let mut a_line = Line2d::new(DVec2::ZERO, DVec2::X); // gp_Lin2d aLine
        let mut a_param = 0.0;
        let mut a_is_valid_segment = the_fexp.segment(a_point, &mut a_line, &mut a_param);
        // TheEdge anEdge — default-constructed (rcad: Option).
        let mut an_edge: Option<Curve2d> = None;
        // OCCT leaves anEdgeOri uninitialized; Forward is the neutral value.
        let mut an_edge_ori = Orientation::Forward;
        let mut a_p_inter = crate::geomalgo::int_res2d::IntersectionPoint::empty();
        let mut a_state = State::Unknown;

        self.nowires = true;

        while a_is_valid_segment {
            self.my_classifier.reset(&a_line, a_param, the_tol);

            // for (theFexp.InitWires(); theFexp.MoreWires(); theFexp.NextWire())
            the_fexp.init_wires();
            while the_fexp.more_wires() {
                self.nowires = false;
                let a_is_w_reject = the_fexp.reject_wire(&a_line, self.my_classifier.parameter());

                if !a_is_w_reject {
                    // Test this wire.
                    // for (theFexp.InitEdges(); theFexp.MoreEdges(); theFexp.NextEdge())
                    the_fexp.init_edges();
                    while the_fexp.more_edges() {
                        let a_is_e_reject =
                            the_fexp.reject_edge(&a_line, self.my_classifier.parameter());

                        if !a_is_e_reject {
                            // Test this edge.
                            // theFexp.CurrentEdge(anEdge, anEdgeOri);
                            let (a_cur_edge, a_cur_ori) = the_fexp.current_edge();
                            an_edge = Some(a_cur_edge.clone());
                            an_edge_ori = a_cur_ori;

                            if an_edge_ori == Orientation::Forward
                                || an_edge_ori == Orientation::Reversed
                            {
                                self.my_classifier.compare(&a_cur_edge, an_edge_ori);
                                let mut a_closest_ind = self.my_classifier.closest_intersection();

                                if a_closest_ind != 0 {
                                    // Save the closest edge.
                                    // auto& anIntersector = theClassifier.Intersector();
                                    let a_nb_pnts =
                                        self.my_classifier.intersector().nb_points();

                                    // theEdge = anEdge;
                                    self.my_edge = an_edge.clone();

                                    if a_closest_ind <= a_nb_pnts {
                                        a_p_inter = self
                                            .my_classifier
                                            .intersector()
                                            .point(a_closest_ind)
                                            .clone();
                                    } else {
                                        a_closest_ind -= a_nb_pnts;

                                        if a_closest_ind & 1 == 1 {
                                            a_p_inter = self
                                                .my_classifier
                                                .intersector()
                                                .segment((a_closest_ind + 1) / 2)
                                                .first_point()
                                                .clone();
                                        } else {
                                            a_p_inter = self
                                                .my_classifier
                                                .intersector()
                                                .segment((a_closest_ind + 1) / 2)
                                                .last_point()
                                                .clone();
                                        }
                                    }

                                    self.my_position =
                                        a_p_inter.transition_of_second().position_on_curve();
                                    self.my_edge_parameter = a_p_inter.param_on_second();
                                }
                                // If we are ON, we stop.
                                a_state = self.my_classifier.state();

                                if a_state == State::On {
                                    return;
                                }
                            }
                        }
                        the_fexp.next_edge();
                    }

                    // If we are out of the wire we stop.
                    a_state = self.my_classifier.state();

                    if a_state == State::Out {
                        return;
                    }
                }
                the_fexp.next_wire();
            }

            if !self.my_classifier.is_head_or_end() && a_state != State::Unknown {
                break;
            }

            // Bad case for classification. Trying to get another segment.
            a_is_valid_segment = the_fexp.other_segment(a_point, &mut a_line, &mut a_param);
        }
    }

    /// OCCT State (cxx L65-68) — the TopClass_FaceClassifier::State body
    /// (TopClass_FaceClassifier.pxx L148-159).
    pub fn state(&self) -> State {
        if self.rejected {
            State::Out
        } else if self.nowires {
            State::In
        } else {
            self.my_classifier.state()
        }
    }

    /// OCCT Rejected (hxx L57) — true when the state was computed by a
    /// rejection (the state is OUT).
    pub fn rejected(&self) -> bool {
        self.rejected
    }

    /// OCCT NoWires (hxx L61) — true if the face contains no wire (the state
    /// is IN).
    pub fn no_wires(&self) -> bool {
        self.nowires
    }

    /// OCCT Edge (cxx L72-75) — the TopClass_FaceClassifier::Edge body
    /// (TopClass_FaceClassifier.pxx L163-168): the edge used to determine
    /// the classification (when the state is ON this is the edge containing
    /// the point).  Raises DomainError after a rejection.
    pub fn edge(&self) -> &Curve2d {
        if self.rejected {
            panic!("TopClass_FaceClassifier::Edge : rejected");
        }
        // OCCT returns the default-constructed adaptor when no edge was
        // recorded; rcad keeps it as None and reports it explicitly.
        self.my_edge
            .as_ref()
            .expect("Geom2dHatch_Classifier::Edge - no edge was recorded")
    }

    /// OCCT EdgeParameter (cxx L79-82) — the TopClass_FaceClassifier::
    /// EdgeParameter body (TopClass_FaceClassifier.pxx L172-176): the
    /// parameter on Edge() used to determine the classification.
    pub fn edge_parameter(&self) -> f64 {
        if self.rejected {
            panic!("TopClass_FaceClassifier::EdgeParameter : rejected");
        }
        self.my_edge_parameter
    }

    /// OCCT Position (hxx L74) — the position of the point on the edge
    /// returned by Edge.
    pub fn position(&self) -> Position {
        self.my_position
    }
}

impl Default for Classifier {
    fn default() -> Self {
        Self::new()
    }
}
