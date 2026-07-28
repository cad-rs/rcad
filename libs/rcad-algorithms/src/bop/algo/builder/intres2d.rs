use glam::DVec2;

// ============================================================================
// IntRes2d_Domain (IntRes2d_Domain.hxx/lxx)
// ============================================================================

/// IntRes2d_Domain — definition of the parameter domain on a 2D curve.
///
/// status bits: bit0=HasFirstPoint, bit1=HasLastPoint, bit2=IsClosed
#[derive(Debug, Clone)]
pub struct IntRes2dDomain {
    status: i32,
    first_param: f64,
    last_param: f64,
    first_tol: f64,
    last_tol: f64,
    first_point: DVec2,
    last_point: DVec2,
    period_first: f64,
    period_last: f64,
}

impl IntRes2dDomain {
    /// Creates an infinite Domain (HasFirstPoint=false, HasLastPoint=false).
    pub fn new() -> Self {
        IntRes2dDomain {
            status: 0,
            first_param: 0.0,
            last_param: 0.0,
            first_tol: 0.0,
            last_tol: 0.0,
            first_point: DVec2::ZERO,
            last_point: DVec2::ZERO,
            period_first: 0.0,
            period_last: 0.0,
        }
    }

    /// Creates a bounded Domain.
    pub fn new_bounded(p1: DVec2, par1: f64, tol1: f64, p2: DVec2, par2: f64, tol2: f64) -> Self {
        let mut d = Self::new();
        d.set_values(p1, par1, tol1, p2, par2, tol2);
        d
    }

    /// OCCT SetValues for a bounded domain.
    pub fn set_values(&mut self, p1: DVec2, par1: f64, tol1: f64, p2: DVec2, par2: f64, tol2: f64) {
        self.status = 3;
        self.period_first = 0.0;
        self.period_last = 0.0;
        self.first_param = par1;
        self.first_point = p1;
        self.first_tol = tol1;
        self.last_param = par2;
        self.last_point = p2;
        self.last_tol = tol2;
    }

    /// Alias for set_values — matches the old rcad name.
    pub fn set_values_bounded(
        &mut self,
        p1: DVec2,
        par1: f64,
        tol1: f64,
        p2: DVec2,
        par2: f64,
        tol2: f64,
    ) {
        self.set_values(p1, par1, tol1, p2, par2, tol2);
    }

    /// OCCT SetValues for an infinite domain.
    pub fn set_values_infinite(&mut self) {
        self.status = 0;
    }

    /// OCCT SetValues for a semi-infinite domain.
    /// If first is true, the point is the first point; otherwise the last.
    pub fn set_values_semi_infinite(&mut self, p: DVec2, par: f64, tol: f64, first: bool) {
        if first {
            self.status = 1;
            self.first_param = par;
            self.first_point = p;
            self.first_tol = tol;
        } else {
            self.status = 2;
            self.last_param = par;
            self.last_point = p;
            self.last_tol = tol;
        }
        self.period_first = 0.0;
        self.period_last = 0.0;
    }

    /// OCCT SetEquivalentParameters — defines a closed domain.
    pub fn set_equivalent_parameters(&mut self, zero: f64, period: f64) {
        assert!((self.status & 3) == 3, "IntRes2dDomain: not bounded");
        self.status |= 4;
        self.period_first = zero;
        self.period_last = period;
    }

    pub fn has_first_point(&self) -> bool {
        (self.status & 1) != 0
    }

    pub fn first_parameter(&self) -> f64 {
        assert!(self.has_first_point(), "IntRes2dDomain: no first point");
        self.first_param
    }

    pub fn first_point(&self) -> DVec2 {
        assert!(self.has_first_point(), "IntRes2dDomain: no first point");
        self.first_point
    }

    pub fn first_tolerance(&self) -> f64 {
        assert!(self.has_first_point(), "IntRes2dDomain: no first point");
        self.first_tol
    }

    pub fn has_last_point(&self) -> bool {
        (self.status & 2) != 0
    }

    pub fn last_parameter(&self) -> f64 {
        assert!(self.has_last_point(), "IntRes2dDomain: no last point");
        self.last_param
    }

    pub fn last_point(&self) -> DVec2 {
        assert!(self.has_last_point(), "IntRes2dDomain: no last point");
        self.last_point
    }

    pub fn last_tolerance(&self) -> f64 {
        assert!(self.has_last_point(), "IntRes2dDomain: no last point");
        self.last_tol
    }

    pub fn is_closed(&self) -> bool {
        (self.status & 4) != 0
    }

    pub fn equivalent_parameters(&self) -> (f64, f64) {
        (self.period_first, self.period_last)
    }
}

// ============================================================================
// IntRes2d_Position / IntRes2d_TypeTrans / IntRes2d_Situation
//   (IntRes2d_Position.hxx, IntRes2d_TypeTrans.hxx, IntRes2d_Situation.hxx)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntRes2dPosition {
    Head,
    Middle,
    End,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntRes2dTypeTrans {
    In,
    Out,
    Touch,
    Undecided,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntRes2dSituation {
    Inside,
    Outside,
    Unknown,
}

// ============================================================================
// IntRes2d_Transition (IntRes2d_Transition.hxx/lxx)
// ============================================================================

/// transition type near an intersection point between two 2D curves.
#[derive(Debug, Clone)]
pub struct IntRes2dTransition {
    tangent: bool,
    posit: IntRes2dPosition,
    typetra: IntRes2dTypeTrans,
    situat: IntRes2dSituation,
    oppos: bool,
}

impl IntRes2dTransition {
    /// Creates an IN or OUT transition.
    pub fn new_in_out(tangent: bool, pos: IntRes2dPosition, typ: IntRes2dTypeTrans) -> Self {
        IntRes2dTransition {
            tangent,
            posit: pos,
            typetra: typ,
            situat: IntRes2dSituation::Unknown,
            oppos: false,
        }
    }

    /// Creates a TOUCH transition.
    pub fn new_touch(
        tangent: bool,
        pos: IntRes2dPosition,
        situ: IntRes2dSituation,
        oppos: bool,
    ) -> Self {
        IntRes2dTransition {
            tangent,
            posit: pos,
            typetra: IntRes2dTypeTrans::Touch,
            situat: situ,
            oppos,
        }
    }

    /// Creates an UNDECIDED transition.
    pub fn new_undecided(pos: IntRes2dPosition) -> Self {
        IntRes2dTransition {
            tangent: true,
            posit: pos,
            typetra: IntRes2dTypeTrans::Undecided,
            situat: IntRes2dSituation::Unknown,
            oppos: false,
        }
    }

    /// OCCT SetValue for IN or OUT.
    pub fn set_value_in_out(
        &mut self,
        tangent: bool,
        pos: IntRes2dPosition,
        typ: IntRes2dTypeTrans,
    ) {
        self.tangent = tangent;
        self.posit = pos;
        self.typetra = typ;
    }

    /// OCCT SetValue for TOUCH.
    pub fn set_value_touch(
        &mut self,
        tangent: bool,
        pos: IntRes2dPosition,
        situ: IntRes2dSituation,
        oppos: bool,
    ) {
        self.tangent = tangent;
        self.posit = pos;
        self.typetra = IntRes2dTypeTrans::Touch;
        self.situat = situ;
        self.oppos = oppos;
    }

    /// OCCT SetValue for UNDECIDED.
    pub fn set_value_undecided(&mut self, pos: IntRes2dPosition) {
        self.posit = pos;
        self.typetra = IntRes2dTypeTrans::Undecided;
    }

    /// OCCT SetPosition.
    pub fn set_position(&mut self, pos: IntRes2dPosition) {
        self.posit = pos;
    }

    /// OCCT PositionOnCurve.
    pub fn position_on_curve(&self) -> IntRes2dPosition {
        self.posit
    }

    /// OCCT TransitionType.
    pub fn transition_type(&self) -> IntRes2dTypeTrans {
        self.typetra
    }

    /// OCCT IsTangent — throws when Undecided.
    pub fn is_tangent(&self) -> bool {
        assert!(
            self.typetra != IntRes2dTypeTrans::Undecided,
            "IntRes2dTransition: IsTangent on Undecided"
        );
        self.tangent
    }

    /// OCCT Situation — throws when not TOUCH.
    pub fn situation(&self) -> IntRes2dSituation {
        assert!(
            self.typetra == IntRes2dTypeTrans::Touch,
            "IntRes2dTransition: Situation on non-Touch"
        );
        self.situat
    }

    /// OCCT IsOpposite — throws when not TOUCH.
    pub fn is_opposite(&self) -> bool {
        assert!(
            self.typetra == IntRes2dTypeTrans::Touch,
            "IntRes2dTransition: IsOpposite on non-Touch"
        );
        self.oppos
    }
}

// ============================================================================
// IntRes2d_IntersectionPoint (IntRes2d_IntersectionPoint.hxx/lxx)
// ============================================================================

/// intersection point between two 2D curves.
#[derive(Debug, Clone)]
pub struct IntRes2dIntersectionPoint {
    pt: DVec2,
    p1: f64,
    p2: f64,
    trans1: IntRes2dTransition,
    trans2: IntRes2dTransition,
}

impl IntRes2dIntersectionPoint {
    /// Creates an IntersectionPoint.
    /// If reversed is true, (Uc1, Trans1) refer to the second curve and
    /// (Uc2, Trans2) refer to the first curve.
    pub fn new(
        p: DVec2,
        u1: f64,
        u2: f64,
        t1: IntRes2dTransition,
        t2: IntRes2dTransition,
        reversed: bool,
    ) -> Self {
        if reversed {
            IntRes2dIntersectionPoint {
                pt: p,
                p1: u2,
                p2: u1,
                trans1: t2,
                trans2: t1,
            }
        } else {
            IntRes2dIntersectionPoint {
                pt: p,
                p1: u1,
                p2: u2,
                trans1: t1,
                trans2: t2,
            }
        }
    }

    /// OCCT SetValues.
    pub fn set_values(
        &mut self,
        p: DVec2,
        u1: f64,
        u2: f64,
        t1: IntRes2dTransition,
        t2: IntRes2dTransition,
        reversed: bool,
    ) {
        if reversed {
            self.pt = p;
            self.p1 = u2;
            self.p2 = u1;
            self.trans1 = t2;
            self.trans2 = t1;
        } else {
            self.pt = p;
            self.p1 = u1;
            self.p2 = u2;
            self.trans1 = t1;
            self.trans2 = t2;
        }
    }

    pub fn value(&self) -> DVec2 {
        self.pt
    }
    pub fn param_on_first(&self) -> f64 {
        self.p1
    }
    pub fn param_on_second(&self) -> f64 {
        self.p2
    }

    /// OCCT TransitionOfFirst.
    pub fn transition_of_first(&self) -> &IntRes2dTransition {
        &self.trans1
    }

    /// OCCT TransitionOfSecond.
    pub fn transition_of_second(&self) -> &IntRes2dTransition {
        &self.trans2
    }
}

// ============================================================================
// IntRes2d_IntersectionSegment (IntRes2d_IntersectionSegment.hxx/lxx)
//   Definition of an intersection curve between two 2D curves.
// ============================================================================

/// segment of intersection between two 2D curves.
#[derive(Debug, Clone)]
pub struct IntRes2dIntersectionSegment {
    oppos: bool,
    first: bool,
    last: bool,
    ptfirst: IntRes2dIntersectionPoint,
    ptlast: IntRes2dIntersectionPoint,
}

impl IntRes2dIntersectionSegment {
    /// Creates an infinite segment (no endpoints).
    pub fn new_infinite(oppos: bool) -> Self {
        let zero_pt = IntRes2dIntersectionPoint {
            pt: DVec2::ZERO,
            p1: 0.0,
            p2: 0.0,
            trans1: IntRes2dTransition::new_undecided(IntRes2dPosition::Middle),
            trans2: IntRes2dTransition::new_undecided(IntRes2dPosition::Middle),
        };
        IntRes2dIntersectionSegment {
            oppos,
            first: false,
            last: false,
            ptfirst: zero_pt.clone(),
            ptlast: zero_pt,
        }
    }

    /// Creates a segment from two endpoints.
    pub fn new_from_points(
        p1: IntRes2dIntersectionPoint,
        p2: IntRes2dIntersectionPoint,
        oppos: bool,
        _reverse_flag: bool,
    ) -> Self {
        IntRes2dIntersectionSegment {
            oppos,
            first: true,
            last: true,
            ptfirst: p1,
            ptlast: p2,
        }
    }

    /// Creates a segment from a single endpoint (semi-infinite).
    pub fn new_from_point(
        p: IntRes2dIntersectionPoint,
        is_first: bool,
        oppos: bool,
        _reverse_flag: bool,
    ) -> Self {
        let zero_pt = IntRes2dIntersectionPoint {
            pt: DVec2::ZERO,
            p1: 0.0,
            p2: 0.0,
            trans1: IntRes2dTransition::new_undecided(IntRes2dPosition::Middle),
            trans2: IntRes2dTransition::new_undecided(IntRes2dPosition::Middle),
        };
        IntRes2dIntersectionSegment {
            oppos,
            first: is_first,
            last: !is_first,
            ptfirst: p,
            ptlast: zero_pt,
        }
    }

    pub fn is_opposite(&self) -> bool {
        self.oppos
    }
    pub fn has_first_point(&self) -> bool {
        self.first
    }
    pub fn has_last_point(&self) -> bool {
        self.last
    }

    pub fn first_point(&self) -> &IntRes2dIntersectionPoint {
        assert!(
            self.has_first_point(),
            "IntRes2dIntersectionSegment: no first point"
        );
        &self.ptfirst
    }

    pub fn last_point(&self) -> &IntRes2dIntersectionPoint {
        assert!(
            self.has_last_point(),
            "IntRes2dIntersectionSegment: no last point"
        );
        &self.ptlast
    }
}
