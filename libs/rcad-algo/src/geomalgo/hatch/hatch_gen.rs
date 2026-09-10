//! OCCT HatchGen package (TKGeomAlgo/HatchGen) — the domain types shared by
//! the Geom2dHatch hatching engine.
//!
//! 1:1 translations:
//!   - HatchGen_IntersectionType.hxx L17-31 — the intersection type enum.
//!   - HatchGen_ErrorStatus.hxx L17-31 — the error status enum.
//!   - HatchGen_IntersectionPoint.hxx L29-94 + .cxx L21-171 — the abstract
//!     base point (C++ inheritance is mapped to the `base` field of
//!     [`PointOnElement`] / [`PointOnHatching`]).
//!   - HatchGen_PointOnElement.hxx L30-69 + .cxx L26-293 + .lxx L22-35.
//!   - HatchGen_PointOnHatching.hxx L31-86 + .cxx L26-237.
//!   - HatchGen_Domain.hxx L27-89 + .cxx L26-95 + .lxx L24-127.

use crate::geomalgo::int_res2d::{
    IntersectionPoint as IntRes2dIntersectionPoint, Position as IntRes2dPosition, Situation,
    TypeTrans,
};
#[cfg(test)]
use crate::geomalgo::int_res2d::Transition;
use rcad_kernel::topods::{Orientation, State};

/// OCCT HatchGen_IntersectionType (hxx L22-28) — intersection type between
/// the hatching and the element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectionType {
    True,
    Touch,
    Tangent,
    Undetermined,
}

/// OCCT HatchGen_ErrorStatus (hxx L21-28) — error status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorStatus {
    NoProblem,
    TrimFailure,
    TransitionFailure,
    IncoherentParity,
    IncompatibleStates,
}

/// OCCT HatchGen_IntersectionPoint (hxx L29-94 + .cxx L21-171) — the base
/// class of HatchGen_PointOnElement and HatchGen_PointOnHatching.  The C++
/// protected members become the fields of this struct; the derived rcad
/// structs embed it as their `base` field.
#[derive(Debug, Clone)]
pub struct IntersectionPoint {
    my_index: i32,
    my_param: f64,
    my_posit: Orientation,
    my_before: State,
    my_after: State,
    my_seg_beg: bool,
    my_seg_end: bool,
}

impl IntersectionPoint {
    /// OCCT HatchGen_IntersectionPoint() (cxx L21-30) — creates an empty
    /// intersection point.
    pub fn new() -> Self {
        IntersectionPoint {
            my_index: 0,
            my_param: f64::MAX, // RealLast()
            my_posit: Orientation::Internal,
            my_before: State::Unknown,
            my_after: State::Unknown,
            my_seg_beg: false,
            my_seg_end: false,
        }
    }

    /// OCCT SetIndex (cxx L37-40).
    pub fn set_index(&mut self, index: i32) {
        self.my_index = index;
    }

    /// OCCT Index (cxx L47-50).
    pub fn index(&self) -> i32 {
        self.my_index
    }

    /// OCCT SetParameter (cxx L57-60).
    pub fn set_parameter(&mut self, parameter: f64) {
        self.my_param = parameter;
    }

    /// OCCT Parameter (cxx L67-70).
    pub fn parameter(&self) -> f64 {
        self.my_param
    }

    /// OCCT SetPosition (cxx L77-80).
    pub fn set_position(&mut self, position: Orientation) {
        self.my_posit = position;
    }

    /// OCCT Position (cxx L87-90).
    pub fn position(&self) -> Orientation {
        self.my_posit
    }

    /// OCCT SetStateBefore (cxx L97-100).
    pub fn set_state_before(&mut self, state: State) {
        self.my_before = state;
    }

    /// OCCT StateBefore (cxx L107-110).
    pub fn state_before(&self) -> State {
        self.my_before
    }

    /// OCCT SetStateAfter (cxx L117-120).
    pub fn set_state_after(&mut self, state: State) {
        self.my_after = state;
    }

    /// OCCT StateAfter (cxx L127-130).
    pub fn state_after(&self) -> State {
        self.my_after
    }

    /// OCCT SetSegmentBeginning (cxx L137-140).
    pub fn set_segment_beginning(&mut self, state: bool) {
        self.my_seg_beg = state;
    }

    /// OCCT SegmentBeginning (cxx L148-151).
    pub fn segment_beginning(&self) -> bool {
        self.my_seg_beg
    }

    /// OCCT SetSegmentEnd (cxx L158-161).
    pub fn set_segment_end(&mut self, state: bool) {
        self.my_seg_end = state;
    }

    /// OCCT SegmentEnd (cxx L168-171).
    pub fn segment_end(&self) -> bool {
        self.my_seg_end
    }
}

impl Default for IntersectionPoint {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT HatchGen_PointOnElement (hxx L30-69 + .cxx L26-293 + .lxx L22-35) —
/// a point on an element.  Inherits HatchGen_IntersectionPoint (rcad: the
/// `base` field) and adds the intersection type.
#[derive(Debug, Clone)]
pub struct PointOnElement {
    /// OCCT HatchGen_IntersectionPoint base class subobject.
    pub base: IntersectionPoint,
    /// OCCT HatchGen_IntersectionType myType.
    my_type: IntersectionType,
}

impl PointOnElement {
    /// OCCT HatchGen_PointOnElement() (cxx L26-29) — creates an empty point
    /// on element.
    pub fn new() -> Self {
        PointOnElement {
            base: IntersectionPoint::new(),
            my_type: IntersectionType::Undetermined,
        }
    }

    /// OCCT HatchGen_PointOnElement(const IntRes2d_IntersectionPoint& Point)
    /// (cxx L33-174) — creates a point from an intersection point.
    pub fn new_intersection(point: &IntRes2dIntersectionPoint) -> Self {
        let trs_h = point.transition_of_first();
        let trs_e = point.transition_of_second();

        let mut r = PointOnElement {
            base: IntersectionPoint::new(),
            my_type: IntersectionType::Undetermined,
        };

        r.base.my_index = 0;
        r.base.my_param = point.param_on_second();

        match trs_e.position_on_curve() {
            IntRes2dPosition::Head => r.base.my_posit = Orientation::Forward,
            IntRes2dPosition::Middle => r.base.my_posit = Orientation::Internal,
            IntRes2dPosition::End => r.base.my_posit = Orientation::Reversed,
        }

        match trs_h.transition_type() {
            TypeTrans::In => {
                r.base.my_before = State::Out;
                r.base.my_after = State::In;
                r.my_type = if r.base.my_posit == Orientation::Internal {
                    IntersectionType::True
                } else {
                    IntersectionType::Touch
                };
            }
            TypeTrans::Out => {
                r.base.my_before = State::In;
                r.base.my_after = State::Out;
                r.my_type = if r.base.my_posit == Orientation::Internal {
                    IntersectionType::True
                } else {
                    IntersectionType::Touch
                };
            }
            // Modified by Sergey KHROMOV - Fri Jan  5 12:07:34 2001 Begin
            TypeTrans::Touch => match trs_h.situation() {
                Situation::Inside => {
                    r.my_type = IntersectionType::Tangent;
                    match r.base.my_posit {
                        Orientation::Forward => {
                            if trs_e.is_opposite() {
                                r.base.my_before = State::In;
                                r.base.my_after = State::Out;
                            } else {
                                r.base.my_before = State::Out;
                                r.base.my_after = State::In;
                            }
                        }
                        Orientation::Internal => {
                            r.base.my_before = State::In;
                            r.base.my_after = State::In;
                        }
                        Orientation::Reversed => {
                            if trs_e.is_opposite() {
                                r.base.my_before = State::Out;
                                r.base.my_after = State::In;
                            } else {
                                r.base.my_before = State::In;
                                r.base.my_after = State::Out;
                            }
                        }
                        Orientation::External => {}
                    }
                }
                Situation::Outside => {
                    r.my_type = IntersectionType::Tangent;
                    match r.base.my_posit {
                        Orientation::Forward => {
                            if trs_e.is_opposite() {
                                r.base.my_before = State::Out;
                                r.base.my_after = State::In;
                            } else {
                                r.base.my_before = State::In;
                                r.base.my_after = State::Out;
                            }
                        }
                        Orientation::Internal => {
                            r.base.my_before = State::Out;
                            r.base.my_after = State::Out;
                        }
                        Orientation::Reversed => {
                            if trs_e.is_opposite() {
                                r.base.my_before = State::In;
                                r.base.my_after = State::Out;
                            } else {
                                r.base.my_before = State::Out;
                                r.base.my_after = State::In;
                            }
                        }
                        Orientation::External => {}
                    }
                }
                Situation::Unknown => {
                    r.base.my_before = State::Unknown;
                    r.base.my_after = State::Unknown;
                    r.my_type = IntersectionType::Tangent;
                }
            },
            // Modified by Sergey KHROMOV - Fri Jan  5 12:07:46 2001 End
            TypeTrans::Undecided => {
                r.base.my_before = State::Unknown;
                r.base.my_after = State::Unknown;
                r.my_type = IntersectionType::Undetermined;
            }
        }

        r.base.my_seg_beg = false;
        r.base.my_seg_end = false;
        r
    }

    // OCCT HatchGen_PointOnElement.lxx L22-35 (inline).

    // -- HatchGen_IntersectionPoint inherited members (the `base` subobject).
    // OCCT exposes them directly on the derived class.

    /// OCCT inherited SetIndex.
    pub fn set_index(&mut self, index: i32) {
        self.base.set_index(index);
    }

    /// OCCT inherited Index.
    pub fn index(&self) -> i32 {
        self.base.index()
    }

    /// OCCT inherited SetParameter.
    pub fn set_parameter(&mut self, parameter: f64) {
        self.base.set_parameter(parameter);
    }

    /// OCCT inherited Parameter.
    pub fn parameter(&self) -> f64 {
        self.base.parameter()
    }

    /// OCCT inherited SetPosition.
    pub fn set_position(&mut self, position: Orientation) {
        self.base.set_position(position);
    }

    /// OCCT inherited Position.
    pub fn position(&self) -> Orientation {
        self.base.position()
    }

    /// OCCT inherited SetStateBefore.
    pub fn set_state_before(&mut self, state: State) {
        self.base.set_state_before(state);
    }

    /// OCCT inherited StateBefore.
    pub fn state_before(&self) -> State {
        self.base.state_before()
    }

    /// OCCT inherited SetStateAfter.
    pub fn set_state_after(&mut self, state: State) {
        self.base.set_state_after(state);
    }

    /// OCCT inherited StateAfter.
    pub fn state_after(&self) -> State {
        self.base.state_after()
    }

    /// OCCT inherited SetSegmentBeginning.
    pub fn set_segment_beginning(&mut self, state: bool) {
        self.base.set_segment_beginning(state);
    }

    /// OCCT inherited SegmentBeginning.
    pub fn segment_beginning(&self) -> bool {
        self.base.segment_beginning()
    }

    /// OCCT inherited SetSegmentEnd.
    pub fn set_segment_end(&mut self, state: bool) {
        self.base.set_segment_end(state);
    }

    /// OCCT inherited SegmentEnd.
    pub fn segment_end(&self) -> bool {
        self.base.segment_end()
    }

    /// OCCT SetIntersectionType (lxx L22-25).
    pub fn set_intersection_type(&mut self, typ: IntersectionType) {
        self.my_type = typ;
    }

    /// OCCT IntersectionType (lxx L32-35).
    pub fn intersection_type(&self) -> IntersectionType {
        self.my_type
    }

    /// OCCT IsIdentical (cxx L181-188) — tests if the point is identical to
    /// an other (all members equal, parameters within Confusion).
    pub fn is_identical(&self, point: &PointOnElement, confusion: f64) -> bool {
        let delta = (self.base.my_param - point.base.my_param).abs();
        (delta <= confusion)
            && (self.base.my_index == point.base.my_index)
            && (self.base.my_posit == point.base.my_posit)
            && (self.my_type == point.my_type)
            && (self.base.my_before == point.base.my_before)
            && (self.base.my_after == point.base.my_after)
            && (self.base.my_seg_beg == point.base.my_seg_beg)
            && (self.base.my_seg_end == point.base.my_seg_end)
    }

    /// OCCT IsDifferent (cxx L195-202) — tests if the point is different from
    /// an other.
    pub fn is_different(&self, point: &PointOnElement, confusion: f64) -> bool {
        let delta = (self.base.my_param - point.base.my_param).abs();
        (delta > confusion)
            || (self.base.my_index != point.base.my_index)
            || (self.base.my_posit != point.base.my_posit)
            || (self.my_type != point.my_type)
            || (self.base.my_before != point.base.my_before)
            || (self.base.my_after != point.base.my_after)
            || (self.base.my_seg_beg != point.base.my_seg_beg)
            || (self.base.my_seg_end != point.base.my_seg_end)
    }
}

impl Default for PointOnElement {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT HatchGen_PointOnHatching (hxx L31-86 + .cxx L26-237) — a point on a
/// hatching carrying the points on the intersected elements.  Inherits
/// HatchGen_IntersectionPoint (rcad: the `base` field).
#[derive(Debug, Clone)]
pub struct PointOnHatching {
    /// OCCT HatchGen_IntersectionPoint base class subobject.
    pub base: IntersectionPoint,
    /// OCCT NCollection_Sequence<HatchGen_PointOnElement> myPoints.
    my_points: Vec<PointOnElement>,
}

impl PointOnHatching {
    /// OCCT HatchGen_PointOnHatching() (cxx L26-28) — creates an empty point.
    pub fn new() -> Self {
        PointOnHatching {
            base: IntersectionPoint::new(),
            my_points: Vec::new(),
        }
    }

    /// OCCT HatchGen_PointOnHatching(const IntRes2d_IntersectionPoint& Point)
    /// (cxx L32-53) — creates a point from an intersection point.
    pub fn new_intersection(point: &IntRes2dIntersectionPoint) -> Self {
        let mut r = PointOnHatching {
            base: IntersectionPoint::new(),
            my_points: Vec::new(),
        };
        r.base.my_index = 0;
        r.base.my_param = point.param_on_first();
        match point.transition_of_first().position_on_curve() {
            IntRes2dPosition::Head => r.base.my_posit = Orientation::Forward,
            IntRes2dPosition::Middle => r.base.my_posit = Orientation::Internal,
            IntRes2dPosition::End => r.base.my_posit = Orientation::Reversed,
        }
        r.base.my_before = State::Unknown;
        r.base.my_after = State::Unknown;
        r.base.my_seg_beg = false;
        r.base.my_seg_end = false;
        r.my_points.clear();
        r
    }

    // -- HatchGen_IntersectionPoint inherited members (the `base` subobject).
    // OCCT exposes them directly on the derived class.

    /// OCCT inherited SetIndex.
    pub fn set_index(&mut self, index: i32) {
        self.base.set_index(index);
    }

    /// OCCT inherited Index.
    pub fn index(&self) -> i32 {
        self.base.index()
    }

    /// OCCT inherited SetParameter.
    pub fn set_parameter(&mut self, parameter: f64) {
        self.base.set_parameter(parameter);
    }

    /// OCCT inherited Parameter.
    pub fn parameter(&self) -> f64 {
        self.base.parameter()
    }

    /// OCCT inherited SetPosition.
    pub fn set_position(&mut self, position: Orientation) {
        self.base.set_position(position);
    }

    /// OCCT inherited Position.
    pub fn position(&self) -> Orientation {
        self.base.position()
    }

    /// OCCT inherited SetStateBefore.
    pub fn set_state_before(&mut self, state: State) {
        self.base.set_state_before(state);
    }

    /// OCCT inherited StateBefore.
    pub fn state_before(&self) -> State {
        self.base.state_before()
    }

    /// OCCT inherited SetStateAfter.
    pub fn set_state_after(&mut self, state: State) {
        self.base.set_state_after(state);
    }

    /// OCCT inherited StateAfter.
    pub fn state_after(&self) -> State {
        self.base.state_after()
    }

    /// OCCT inherited SetSegmentBeginning.
    pub fn set_segment_beginning(&mut self, state: bool) {
        self.base.set_segment_beginning(state);
    }

    /// OCCT inherited SegmentBeginning.
    pub fn segment_beginning(&self) -> bool {
        self.base.segment_beginning()
    }

    /// OCCT inherited SetSegmentEnd.
    pub fn set_segment_end(&mut self, state: bool) {
        self.base.set_segment_end(state);
    }

    /// OCCT inherited SegmentEnd.
    pub fn segment_end(&self) -> bool {
        self.base.segment_end()
    }

    /// OCCT AddPoint (cxx L60-74) — adds a point on element to the point.
    pub fn add_point(&mut self, point: &PointOnElement, confusion: f64) {
        let nb_pnt = self.my_points.len() as i32;
        // for (IPnt = 1; IPnt <= NbPnt && myPoints(IPnt).IsDifferent(Point,
        // Confusion); IPnt++) { ; }
        let mut ipnt = 1i32;
        while ipnt <= nb_pnt
            && self.my_points[(ipnt - 1) as usize].is_different(point, confusion)
        {
            ipnt += 1;
        }
        if ipnt > nb_pnt {
            self.my_points.push(point.clone());
        }
    }

    /// OCCT NbPoints (cxx L82-85) — the number of elements intersecting the
    /// hatching at this point.
    pub fn nb_points(&self) -> usize {
        self.my_points.len()
    }

    /// OCCT Point (cxx L92-95) — returns the Index-th point on element of the
    /// point (1-based; OCCT raises OutOfRange if Index > NbPoints).
    pub fn point(&self, index: usize) -> &PointOnElement {
        &self.my_points[index - 1]
    }

    /// OCCT RemPoint (cxx L102-105) — removes the Index-th point on element
    /// of the point (1-based).
    pub fn rem_point(&mut self, index: usize) {
        self.my_points.remove(index - 1);
    }

    /// OCCT ClrPoints (cxx L112-115) — removes all the points on element of
    /// the point.
    pub fn clr_points(&mut self) {
        self.my_points.clear();
    }

    /// OCCT IsLower (cxx L122-126) — a point on hatching P1 is lower than P2
    /// if P2.myParam - P1.myParam > Confusion.
    pub fn is_lower(&self, point: &PointOnHatching, confusion: f64) -> bool {
        point.base.my_param - self.base.my_param > confusion
    }

    /// OCCT IsEqual (cxx L133-137) — |P2.myParam - P1.myParam| <= Confusion.
    pub fn is_equal(&self, point: &PointOnHatching, confusion: f64) -> bool {
        (point.base.my_param - self.base.my_param).abs() <= confusion
    }

    /// OCCT IsGreater (cxx L144-148) — P1.myParam - P2.myParam > Confusion.
    pub fn is_greater(&self, point: &PointOnHatching, confusion: f64) -> bool {
        self.base.my_param - point.base.my_param > confusion
    }
}

impl Default for PointOnHatching {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT HatchGen_Domain (hxx L27-89 + .cxx L26-95 + .lxx L24-127) — the
/// trimmed domain of a hatching curve, bounded by two points on the hatching
/// (either side may be infinite).
#[derive(Debug, Clone)]
pub struct Domain {
    my_has_first_point: bool,
    my_first_point: PointOnHatching,
    my_has_second_point: bool,
    my_second_point: PointOnHatching,
}

impl Domain {
    /// OCCT HatchGen_Domain() (cxx L26-30) — creates an infinite domain.
    pub fn new() -> Self {
        Domain {
            my_has_first_point: false,
            my_first_point: PointOnHatching::new(),
            my_has_second_point: false,
            my_second_point: PointOnHatching::new(),
        }
    }

    /// OCCT HatchGen_Domain(P1, P2) (cxx L34-41) — creates a domain for the
    /// curve associated to a hatching.
    pub fn new_bounded(p1: &PointOnHatching, p2: &PointOnHatching) -> Self {
        Domain {
            my_has_first_point: true,
            my_first_point: p1.clone(),
            my_has_second_point: true,
            my_second_point: p2.clone(),
        }
    }

    /// OCCT HatchGen_Domain(P, First) (cxx L45-59) — creates a semi-infinite
    /// domain. The `first` flag means that the given point is the first one.
    pub fn new_semi(p: &PointOnHatching, first: bool) -> Self {
        let mut r = Domain {
            my_has_first_point: false,
            my_first_point: PointOnHatching::new(),
            my_has_second_point: false,
            my_second_point: PointOnHatching::new(),
        };
        if first {
            r.my_has_first_point = true;
            r.my_has_second_point = false;
            r.my_first_point = p.clone();
        } else {
            r.my_has_first_point = false;
            r.my_has_second_point = true;
            r.my_second_point = p.clone();
        }
        r
    }

    // OCCT HatchGen_Domain.lxx L24-127 (inline methods).

    /// OCCT SetPoints(P1, P2) (lxx L24-31).
    pub fn set_points(&mut self, p1: &PointOnHatching, p2: &PointOnHatching) {
        self.my_has_first_point = true;
        self.my_first_point = p1.clone();
        self.my_has_second_point = true;
        self.my_second_point = p2.clone();
    }

    /// OCCT SetPoints() (lxx L39-43) — sets both points at the infinite.
    pub fn set_points_infinite(&mut self) {
        self.my_has_first_point = false;
        self.my_has_second_point = false;
    }

    /// OCCT SetFirstPoint(P) (lxx L50-54).
    pub fn set_first_point(&mut self, p: &PointOnHatching) {
        self.my_has_first_point = true;
        self.my_first_point = p.clone();
    }

    /// OCCT SetFirstPoint() (lxx L61-64) — sets the first point at the
    /// infinite.
    pub fn set_first_point_infinite(&mut self) {
        self.my_has_first_point = false;
    }

    /// OCCT SetSecondPoint(P) (lxx L71-75).
    pub fn set_second_point(&mut self, p: &PointOnHatching) {
        self.my_has_second_point = true;
        self.my_second_point = p.clone();
    }

    /// OCCT SetSecondPoint() (lxx L82-85) — sets the second point at the
    /// infinite.
    pub fn set_second_point_infinite(&mut self) {
        self.my_has_second_point = false;
    }

    /// OCCT HasFirstPoint (lxx L92-95).
    pub fn has_first_point(&self) -> bool {
        self.my_has_first_point
    }

    /// OCCT FirstPoint (lxx L102-106) — raises DomainError when the domain
    /// has no first point.
    pub fn first_point(&self) -> &PointOnHatching {
        if !self.my_has_first_point {
            panic!("HatchGen_Domain::FirstPoint");
        }
        &self.my_first_point
    }

    /// OCCT HasSecondPoint (lxx L113-116).
    pub fn has_second_point(&self) -> bool {
        self.my_has_second_point
    }

    /// OCCT SecondPoint (lxx L123-127) — raises DomainError when the domain
    /// has no second point.
    pub fn second_point(&self) -> &PointOnHatching {
        if !self.my_has_second_point {
            panic!("HatchGen_Domain::SecondPoint");
        }
        &self.my_second_point
    }
}

impl Default for Domain {
    fn default() -> Self {
        Self::new()
    }
}

/// Helper used by the tests below: an IntRes2d point with one IN/OUT
/// transition on the first curve and a plain Middle position on the second.
#[cfg(test)]
pub(crate) fn make_intres2d_point(
    p: glam::DVec2,
    par1: f64,
    par2: f64,
    trs1: Transition,
    trs2: Transition,
) -> IntRes2dIntersectionPoint {
    IntRes2dIntersectionPoint::new(p, par1, par2, trs1, trs2, false)
}
