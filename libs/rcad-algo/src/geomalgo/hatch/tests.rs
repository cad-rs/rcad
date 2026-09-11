//! Anchor tests for the Geom2dHatch / HatchGen leaf layer.  Each leaf gets
//! at least one analytically assertable OCCT anchor (line/circle elements,
//! IN/OUT/ON classification, line x line intersections).

use super::classifier::Classifier;
use super::element::HatchElement;
use super::elements::HatchElements;
use super::fclass2d::FClass2dOfClassifier;
use super::hatch_gen::{
    Domain, ErrorStatus, IntersectionType, PointOnElement, PointOnHatching,
};
use super::hatching::Hatching;
use super::intersector::HatchIntersector;
use crate::geomalgo::int_res2d::{Position, Situation, Transition, TypeTrans};
use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Line2d, TrimmedCurve2};
use rcad_kernel::topods::{Orientation, State};

/// An IntRes2d_IntersectionPoint at the given parameters/transitions (the
/// HatchGen test helper).
fn make_point(
    p: DVec2,
    par1: f64,
    par2: f64,
    trs1: Transition,
    trs2: Transition,
) -> crate::geomalgo::int_res2d::IntersectionPoint {
    crate::geomalgo::int_res2d::IntersectionPoint::new(p, par1, par2, trs1, trs2, false)
}

/// A bounded element on the segment [p1, p2] (the rcad analogue of a
/// Geom2dAdaptor_Curve over a trimmed Geom2d_Line).
fn seg_element(p1: DVec2, p2: DVec2, orientation: Orientation) -> HatchElement {
    let dir = (p2 - p1).normalize();
    HatchElement::new(
        Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Line(Line2d::new(p1, dir))),
            t_min: 0.0,
            t_max: 1.0,
        }),
        orientation,
    )
}

/// The unit square boundary, CCW, all edges FORWARD (material on the left).
fn bind_square(elements: &mut HatchElements) {
    elements.bind(1, seg_element(DVec2::new(0., 0.), DVec2::new(1., 0.), Orientation::Forward));
    elements.bind(2, seg_element(DVec2::new(1., 0.), DVec2::new(1., 1.), Orientation::Forward));
    elements.bind(3, seg_element(DVec2::new(1., 1.), DVec2::new(0., 1.), Orientation::Forward));
    elements.bind(4, seg_element(DVec2::new(0., 1.), DVec2::new(0., 0.), Orientation::Forward));
}

// ---------------------------------------------------------------------------
// HatchGen_IntersectionPoint / PointOnElement — the state machine anchor.
// ---------------------------------------------------------------------------

#[test]
fn point_on_element_state_machine() {
    // A crossing point: the line goes IN (Middle/Middle) -> TRUE, OUT->IN,
    // position INTERNAL.
    let trs_h = Transition::in_out(false, Position::Middle, TypeTrans::In);
    let trs_e = Transition::in_out(false, Position::Middle, TypeTrans::Out);
    let p = make_point(
        DVec2::new(0.5, 0.0),
        0.629,
        0.877,
        trs_h,
        trs_e,
    );
    let pOE = PointOnElement::new_intersection(&p);
    assert_eq!(pOE.index(), 0);
    assert!((pOE.parameter() - 0.877).abs() < 1e-12); // ParamOnSecond
    assert_eq!(pOE.position(), Orientation::Internal); // Middle
    assert_eq!(pOE.intersection_type(), IntersectionType::True);
    assert_eq!(pOE.state_before(), State::Out);
    assert_eq!(pOE.state_after(), State::In);
    assert!(!pOE.segment_beginning());
    assert!(!pOE.segment_end());

    // A tangent touch at the head of the edge (Inside, not opposite):
    // TANGENT with IN/OUT for FORWARD position.
    let trs_h_touch = Transition::touch(false, Position::Middle, Situation::Inside, false);
    let trs_e_head = Transition::in_out(false, Position::Head, TypeTrans::In);
    let p2 = make_point(
        DVec2::ZERO,
        0.0,
        0.0,
        trs_h_touch,
        trs_e_head,
    );
    let pOE2 = PointOnElement::new_intersection(&p2);
    assert_eq!(pOE2.position(), Orientation::Forward); // Head
    assert_eq!(pOE2.intersection_type(), IntersectionType::Tangent);
    assert_eq!(pOE2.state_before(), State::Out);
    assert_eq!(pOE2.state_after(), State::In);

    // An UNDECIDED transition -> UNDETERMINED with UNKNOWN states.
    let trs_undecided = Transition::undecided(Position::Middle);
    let p3 = make_point(
        DVec2::ZERO,
        1.0,
        2.0,
        trs_undecided,
        trs_undecided,
    );
    let pOE3 = PointOnElement::new_intersection(&p3);
    assert_eq!(pOE3.intersection_type(), IntersectionType::Undetermined);
    assert_eq!(pOE3.state_before(), State::Unknown);
    assert_eq!(pOE3.state_after(), State::Unknown);

    // IsIdentical / IsDifferent.
    assert!(pOE.is_identical(&pOE.clone(), 1e-7));
    assert!(!pOE.is_different(&pOE.clone(), 1e-7));
    let mut pOE4 = pOE.clone();
    pOE4.set_parameter(pOE.parameter() + 1e-3);
    assert!(pOE.is_different(&pOE4, 1e-7));
    assert!(!pOE.is_identical(&pOE4, 1e-7));
}

// ---------------------------------------------------------------------------
// HatchGen_PointOnHatching + HatchGen_Domain.
// ---------------------------------------------------------------------------

#[test]
fn point_on_hatching_and_domain() {
    // PointOnHatching takes ParamOnFirst.
    let trs_h = Transition::in_out(false, Position::Middle, TypeTrans::In);
    let trs_e = Transition::in_out(false, Position::Middle, TypeTrans::Out);
    let p = make_point(DVec2::ZERO, 0.25, 0.75, trs_h, trs_e);
    let pOH = PointOnHatching::new_intersection(&p);
    assert!((pOH.parameter() - 0.25).abs() < 1e-12);

    // AddPoint deduplicates identical points within Confusion.
    let pOE1 = PointOnElement::new_intersection(&p);
    let mut pOH2 = PointOnHatching::new();
    pOH2.set_parameter(0.25);
    pOH2.add_point(&pOE1, 1e-7);
    pOH2.add_point(&pOE1, 1e-7);
    assert_eq!(pOH2.nb_points(), 1);

    // A different point is appended.
    let q = make_point(
        DVec2::ZERO,
        0.25,
        0.5,
        Transition::in_out(false, Position::Middle, TypeTrans::Out),
        trs_e,
    );
    let pOE2 = PointOnElement::new_intersection(&q);
    pOH2.add_point(&pOE2, 1e-7);
    assert_eq!(pOH2.nb_points(), 2);

    // Lower / equal / greater with Confusion = 1e-3.
    let mut other = PointOnHatching::new();
    other.set_parameter(0.25);
    assert!(pOH2.is_equal(&other, 1e-3));
    let mut greater = PointOnHatching::new();
    greater.set_parameter(0.25 - 1.0);
    assert!(pOH2.is_greater(&greater, 1e-3)); // 0.25 - (0.25-1) = 1 > 1e-3
    let mut lower = PointOnHatching::new();
    lower.set_parameter(0.25 + 1.0);
    assert!(pOH2.is_lower(&lower, 1e-3));

    // RemPoint / ClrPoints.
    pOH2.rem_point(2);
    assert_eq!(pOH2.nb_points(), 1);
    pOH2.clr_points();
    assert_eq!(pOH2.nb_points(), 0);

    // Domain state machine: bounded / semi / infinite.
    let mut d = Domain::new_bounded(&pOH, &pOH2);
    assert!(d.has_first_point() && d.has_second_point());
    d.set_second_point_infinite();
    assert!(!d.has_second_point());
    d.set_points_infinite();
    assert!(!d.has_first_point());
    // FirstPoint raises DomainError without a first point.
    let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = d.first_point();
    }));
    assert!(panicked.is_err());
    d.set_first_point(&pOH);
    assert_eq!(d.first_point().parameter(), 0.25);
    let semi = Domain::new_semi(&pOH, false); // the given point is the second
    assert!(!semi.has_first_point() && semi.has_second_point());
}

// ---------------------------------------------------------------------------
// Geom2dHatch_Element / Elements — add/remove/lookup + traversal.
// ---------------------------------------------------------------------------

#[test]
fn elements_bind_find_traverse() {
    let mut elements = HatchElements::new();
    assert!(!elements.is_bound(1));
    elements.bind(1, seg_element(DVec2::ZERO, DVec2::X, Orientation::Forward));
    elements.bind(2, seg_element(DVec2::X, DVec2::new(1., 1.), Orientation::Reversed));
    assert!(elements.is_bound(1) && elements.is_bound(2));

    // Bind over an existing key overrides and returns false.
    assert!(!elements.bind(1, seg_element(DVec2::ZERO, DVec2::X, Orientation::Reversed)));
    assert_eq!(elements.find(1).orientation(), Orientation::Reversed);

    // ChangeFind mutates the bound element.
    elements.change_find(2).set_orientation(Orientation::Forward);
    assert_eq!(elements.find(2).orientation(), Orientation::Forward);

    // Edge traversal over the single wire.
    elements.init_wires();
    assert!(elements.more_wires());
    elements.next_wire();
    assert!(!elements.more_wires());

    elements.init_edges();
    let mut n = 0;
    let mut orientations = Vec::new();
    while elements.more_edges() {
        let (_curve, ori) = elements.current_edge();
        orientations.push(ori);
        n += 1;
        elements.next_edge();
    }
    assert_eq!(n, 2);
    assert_eq!(orientations, vec![Orientation::Reversed, Orientation::Forward]);

    // UnBind / Clear.
    assert!(elements.un_bind(1));
    assert!(!elements.is_bound(1));
    assert!(!elements.un_bind(1)); // nothing left to unbind
    elements.clear();
    assert!(!elements.is_bound(2));

    // The OCCT copy constructor is "magic": an empty map with reset state
    // (OCC12627) — cloning does not copy the source elements.
    let mut src = HatchElements::new();
    src.bind(1, seg_element(DVec2::ZERO, DVec2::X, Orientation::Forward));
    let copy = src.clone();
    assert!(!copy.is_bound(1));
}

#[test]
fn elements_segment_probing_on_square() {
    let mut elements = HatchElements::new();
    bind_square(&mut elements);

    // From the center: the probing must produce a line through P crossing
    // the boundary, with the parameter = distance to the boundary point.
    let mut line = Line2d::new(DVec2::ZERO, DVec2::X);
    let mut par = 0.0;
    assert!(elements.segment(DVec2::new(0.5, 0.5), &mut line, &mut par));
    assert!(par > 0.1 && par < 2.0);
    // The produced line passes through P (its location) — sanity anchor.
    assert!((line.origin - DVec2::new(0.5, 0.5)).length() < 1e-12);

    // From a far outside point: probing still finds a segment.
    let mut elements2 = HatchElements::new();
    bind_square(&mut elements2);
    let mut par2 = 0.0;
    assert!(elements2.other_segment(DVec2::new(3.0, 3.0), &mut line, &mut par2));
    assert!(par2 > 0.1);

    // An empty element set yields no segment (Par = RealLast).
    let mut empty = HatchElements::new();
    let mut par3 = 0.0;
    assert!(!empty.segment(DVec2::new(0.5, 0.5), &mut line, &mut par3));
    assert_eq!(par3, f64::MAX);
}

// ---------------------------------------------------------------------------
// Geom2dHatch_Intersector — line x line intersections.
// ---------------------------------------------------------------------------

#[test]
fn intersector_perform_line_line() {
    let mut inter = HatchIntersector::with_tolerances(1.0e-7, 1.0e-7);
    assert_eq!(inter.confusion_tolerance(), 1.0e-7);
    assert_eq!(inter.tangency_tolerance(), 1.0e-7);

    // The hatching segment: from (0,0) along +x, parameter range [0, 2].
    let l = Line2d::new(DVec2::ZERO, DVec2::X);

    // A crossing vertical edge: (0.5,-1) -> (0.5,1); the trimmed curve keeps
    // the basis parameterization, so the crossing point sits at t = 0.
    let cross = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(Line2d::new(DVec2::new(0.5, 0.0), DVec2::Y))),
        t_min: -1.0,
        t_max: 1.0,
    });
    inter.perform(&l, 2.0, 1.0e-7, &cross);
    assert!(inter.is_done());
    assert_eq!(inter.nb_points(), 1);
    assert_eq!(inter.nb_segments(), 0);
    let p = inter.point(1);
    assert!((p.param_on_first() - 0.5).abs() < 1e-9);
    assert!(p.param_on_second().abs() < 1e-9);

    // A parallel edge: no intersection at all.
    let parallel = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(Line2d::new(DVec2::new(-1.0, 1.0), DVec2::X))),
        t_min: -1.0,
        t_max: 1.0,
    });
    let mut inter2 = HatchIntersector::with_tolerances(1.0e-7, 1.0e-7);
    inter2.perform(&l, 2.0, 1.0e-7, &parallel);
    assert_eq!(inter2.nb_points(), 0);
    assert_eq!(inter2.nb_segments(), 0);
}

// ---------------------------------------------------------------------------
// Geom2dHatch_FClass2dOfClassifier — a single-edge classification anchor.
// ---------------------------------------------------------------------------

#[test]
fn fclass2d_compare_crossing_state() {
    let mut cl = FClass2dOfClassifier::new();
    assert_eq!(cl.state(), State::Unknown);
    assert_eq!(cl.closest_intersection(), 0);

    // The probing segment: from the center of the square leftwards.
    let line = Line2d::new(DVec2::new(0.5, 0.5), DVec2::new(-1.0, 0.0));
    cl.reset(&line, 5.0, 1.0e-7);

    // The left edge (FORWARD, material on +x): the ray EXITS the material
    // through it, so the point state is IN.
    let left = seg_element(DVec2::new(0., 1.), DVec2::new(0., 0.), Orientation::Forward);
    cl.compare(left.curve(), Orientation::Forward);
    assert_eq!(cl.closest_intersection(), 1);
    assert_eq!(cl.state(), State::In);
    assert!(!cl.is_head_or_end());

    // A parallel edge produces no intersection: the state machine keeps the
    // previous state but drops the closest intersection.
    let bottom = seg_element(DVec2::ZERO, DVec2::X, Orientation::Forward);
    cl.compare(bottom.curve(), Orientation::Forward);
    assert_eq!(cl.closest_intersection(), 0);
    assert_eq!(cl.state(), State::In);
}

// ---------------------------------------------------------------------------
// Geom2dHatch_Classifier — IN / OUT / ON on a simple closed domain.
// ---------------------------------------------------------------------------

#[test]
fn classifier_square_in_out_on() {
    // Inside point.
    let mut elements = HatchElements::new();
    bind_square(&mut elements);
    let mut cl = Classifier::new();
    cl.perform(&mut elements, DVec2::new(0.5, 0.5), 1.0e-7);
    assert!(!cl.rejected());
    assert!(!cl.no_wires());
    assert_eq!(cl.state(), State::In);

    // Outside point.
    let mut elements2 = HatchElements::new();
    bind_square(&mut elements2);
    let mut cl2 = Classifier::new();
    cl2.perform(&mut elements2, DVec2::new(2.0, 2.0), 1.0e-7);
    assert_eq!(cl2.state(), State::Out);

    // Point on the bottom edge: ON with the recorded edge/parameter.
    let mut elements3 = HatchElements::new();
    bind_square(&mut elements3);
    let mut cl3 = Classifier::new();
    cl3.perform(&mut elements3, DVec2::new(0.5, 0.0), 1.0e-7);
    assert_eq!(cl3.state(), State::On);
    assert_eq!(cl3.position(), Position::Middle);
    assert!((cl3.edge_parameter() - 0.5).abs() < 1e-6);
}

// ---------------------------------------------------------------------------
// Geom2dHatch_Hatching — points / domains bookkeeping.
// ---------------------------------------------------------------------------

#[test]
fn hatching_points_and_domains() {
    let mut h = Hatching::new(Curve2d::Line(Line2d::new(DVec2::ZERO, DVec2::X)));
    assert_eq!(h.status(), ErrorStatus::NoProblem);
    assert!(!h.is_done());
    assert!(!h.trim_done() && !h.trim_failed());

    // ClassificationPoint of an infinite-range line: a = -inf, b = +inf -> 0.
    let p0 = h.classification_point();
    assert_eq!(p0, DVec2::ZERO);

    // TrimFailed propagates into the status.
    h.set_trim_failed(true);
    assert!(h.trim_failed());
    assert_eq!(h.status(), ErrorStatus::TrimFailure);
    h.set_status(ErrorStatus::NoProblem);
    h.set_trim_failed(false);

    // Points.
    let trs = Transition::in_out(false, Position::Middle, TypeTrans::In);
    let ip = make_point(
        DVec2::ZERO,
        0.5,
        0.5,
        trs,
        Transition::in_out(false, Position::Middle, TypeTrans::Out),
    );
    let pOH = PointOnHatching::new_intersection(&ip);
    h.add_point(&pOH, 1e-7);
    h.add_point(&pOH, 1e-7);
    assert_eq!(h.nb_points(), 1);
    assert!((h.point(1).parameter() - 0.5).abs() < 1e-12);

    // Domains: IsDone flag invalidation on point mutations.
    h.set_is_done(true);
    let d = Domain::new();
    h.add_domain(&d);
    assert_eq!(h.nb_domains(), 1);
    // add_point clears the domains when IsDone.
    h.add_point(&pOH, 1e-7);
    assert_eq!(h.nb_domains(), 0);
    assert!(!h.is_done());

    h.set_is_done(true);
    h.rem_point(1);
    assert_eq!(h.nb_points(), 0);
    h.clr_points();
    assert!(!h.trim_done() && !h.trim_failed());

    // ChangeCurve / curve accessors.  A plain line spans the infinite domain
    // (Geom2d_Line.cxx L142-151: -/+ Precision::Infinite()).
    h.change_curve();
    assert_eq!(
        h.curve().default_domain(),
        [
            -rcad_kernel::core::precision::INFINITE_VALUE,
            rcad_kernel::core::precision::INFINITE_VALUE
        ]
    );
}
