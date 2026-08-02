// OCCT IntRes2d package (IntRes2d_Transition.hxx/.lxx/.cxx,
// IntRes2d_Domain.hxx/.lxx/.cxx, IntRes2d_IntersectionPoint.hxx/.lxx/.cxx,
// IntRes2d_IntersectionSegment.hxx/.lxx/.cxx)
//
// Core data types for 2D curve intersection (used by BRepClass_Intersector
// and Geom2dInt_GInter in the FClass2d face-classifier fallback).

use glam::DVec2;

/// OCCT IntRes2d_Position — where on the curve an intersection occurs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Position {
    Head,
    Middle,
    End,
}

/// OCCT IntRes2d_TypeTrans — type of transition near an intersection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeTrans {
    In,
    Out,
    Touch,
    Undecided,
}

/// OCCT IntRes2d_Situation — TOUCH-transition situation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Situation {
    Inside,
    Outside,
    Unknown,
}

/// OCCT IntRes2d_Transition — transition of one curve near an intersection
/// with another.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transition {
    tangent: bool,
    posit: Position,
    typetra: TypeTrans,
    situat: Situation,
    oppos: bool,
}

impl Transition {
    /// Empty constructor.
    pub fn empty() -> Self {
        Transition {
            tangent: false,
            posit: Position::Middle,
            typetra: TypeTrans::Undecided,
            situat: Situation::Unknown,
            oppos: false,
        }
    }

    /// OCCT IN/OUT transition: IntRes2d_Transition(Tangent, Pos, Type).
    pub fn in_out(tangent: bool, pos: Position, typetra: TypeTrans) -> Self {
        Transition {
            tangent,
            posit: pos,
            typetra,
            situat: Situation::Unknown,
            oppos: false,
        }
    }

    /// OCCT TOUCH transition: IntRes2d_Transition(Tangent, Pos, Situ, Oppos).
    pub fn touch(tangent: bool, pos: Position, situat: Situation, oppos: bool) -> Self {
        Transition {
            tangent,
            posit: pos,
            typetra: TypeTrans::Touch,
            situat,
            oppos,
        }
    }

    /// OCCT UNDECIDED transition: IntRes2d_Transition(Pos).
    pub fn undecided(pos: Position) -> Self {
        Transition {
            tangent: true,
            posit: pos,
            typetra: TypeTrans::Undecided,
            situat: Situation::Unknown,
            oppos: false,
        }
    }

    pub fn set_position(&mut self, pos: Position) {
        self.posit = pos;
    }

    /// OCCT PositionOnCurve().
    pub fn position_on_curve(&self) -> Position {
        self.posit
    }

    /// OCCT TransitionType().
    pub fn transition_type(&self) -> TypeTrans {
        self.typetra
    }

    /// OCCT IsTangent().
    pub fn is_tangent(&self) -> bool {
        self.tangent
    }

    /// OCCT Situation().
    pub fn situation(&self) -> Situation {
        self.situat
    }

    /// OCCT IsOpposite().
    pub fn is_opposite(&self) -> bool {
        self.oppos
    }
}

/// OCCT IntRes2d_Domain — domain of parameter on a 2D curve.
///
/// status bit flags (IntRes2d_Domain.lxx L19-23):
/// - has first point ↔ status & 1
/// - has last point  ↔ status & 2
/// - closed          ↔ status & 4
#[derive(Debug, Clone)]
pub struct Domain {
    status: i32,
    first_param: f64,
    last_param: f64,
    first_tol: f64,
    last_tol: f64,
    first_point: DVec2,
    last_point: DVec2,
    periodfirst: f64,
    periodlast: f64,
}

// OCCT IntRes2d_Domain.cxx L28-31: LimitInfinite
fn limit_infinite(val: f64) -> f64 {
    const INF_VAL: f64 = 1.0e100;
    if val.abs() > INF_VAL {
        if val > 0.0 {
            INF_VAL
        } else {
            -INF_VAL
        }
    } else {
        val
    }
}

impl Domain {
    /// Creates an infinite Domain (no first/last point).
    pub fn infinite() -> Self {
        Domain {
            status: 0,
            first_param: 0.0,
            last_param: 0.0,
            first_tol: 0.0,
            last_tol: 0.0,
            first_point: DVec2::ZERO,
            last_point: DVec2::ZERO,
            periodfirst: 0.0,
            periodlast: 0.0,
        }
    }

    /// Creates a bounded Domain.
    pub fn bounded(
        pnt1: DVec2,
        par1: f64,
        tol1: f64,
        pnt2: DVec2,
        par2: f64,
        tol2: f64,
    ) -> Self {
        let mut d = Domain::infinite();
        d.set_values_bounded(pnt1, par1, tol1, pnt2, par2, tol2);
        d
    }

    /// Creates a semi-infinite Domain. If `first` is true, the point is the
    /// first point; otherwise the last.
    pub fn semi(pnt: DVec2, par: f64, tol: f64, first: bool) -> Self {
        let mut d = Domain::infinite();
        d.set_values_semi(pnt, par, tol, first);
        d
    }

    /// OCCT SetValues() — infinite domain.
    pub fn set_values(&mut self) {
        self.status = 0;
        self.periodfirst = 0.0;
        self.periodlast = 0.0;
    }

    /// OCCT SetValues(Pnt1, Par1, Tol1, Pnt2, Par2, Tol2) — bounded domain.
    pub fn set_values_bounded(
        &mut self,
        pnt1: DVec2,
        par1: f64,
        tol1: f64,
        pnt2: DVec2,
        par2: f64,
        tol2: f64,
    ) {
        self.status = 3;
        self.periodfirst = 0.0;
        self.periodlast = 0.0;
        self.first_param = limit_infinite(par1);
        self.first_point = DVec2::new(limit_infinite(pnt1.x), limit_infinite(pnt1.y));
        self.first_tol = tol1;
        self.last_param = limit_infinite(par2);
        self.last_point = DVec2::new(limit_infinite(pnt2.x), limit_infinite(pnt2.y));
        self.last_tol = tol2;
    }

    /// OCCT SetValues(Pnt, Par, Tol, First) — semi-infinite domain.
    pub fn set_values_semi(&mut self, pnt: DVec2, par: f64, tol: f64, first: bool) {
        self.periodfirst = 0.0;
        self.periodlast = 0.0;
        if first {
            self.status = 1;
            self.first_param = limit_infinite(par);
            self.first_point = DVec2::new(limit_infinite(pnt.x), limit_infinite(pnt.y));
            self.first_tol = tol;
        } else {
            self.status = 2;
            self.last_param = limit_infinite(par);
            self.last_point = DVec2::new(limit_infinite(pnt.x), limit_infinite(pnt.y));
            self.last_tol = tol;
        }
    }

    /// OCCT SetEquivalentParameters(zero, period) — mark the domain closed.
    pub fn set_equivalent_parameters(&mut self, p_first: f64, p_last: f64) {
        self.status |= 4;
        self.periodfirst = p_first;
        self.periodlast = p_last;
    }

    pub fn has_first_point(&self) -> bool {
        (self.status & 1) != 0
    }

    pub fn first_parameter(&self) -> f64 {
        self.first_param
    }

    pub fn first_point(&self) -> DVec2 {
        self.first_point
    }

    pub fn first_tolerance(&self) -> f64 {
        self.first_tol
    }

    pub fn has_last_point(&self) -> bool {
        (self.status & 2) != 0
    }

    pub fn last_parameter(&self) -> f64 {
        self.last_param
    }

    pub fn last_point(&self) -> DVec2 {
        self.last_point
    }

    pub fn last_tolerance(&self) -> f64 {
        self.last_tol
    }

    pub fn is_closed(&self) -> bool {
        (self.status & 4) != 0
    }

    pub fn equivalent_parameters(&self) -> (f64, f64) {
        (self.periodfirst, self.periodlast)
    }
}

/// OCCT IntRes2d_IntersectionPoint — an intersection point between two 2D
/// curves, with the parameter on each curve and the transition of each.
#[derive(Debug, Clone)]
pub struct IntersectionPoint {
    pt: DVec2,
    p1: f64,
    p2: f64,
    trans1: Transition,
    trans2: Transition,
}

impl IntersectionPoint {
    pub fn empty() -> Self {
        IntersectionPoint {
            pt: DVec2::ZERO,
            p1: f64::MAX, // RealLast
            p2: f64::MAX,
            trans1: Transition::empty(),
            trans2: Transition::empty(),
        }
    }

    /// OCCT IntRes2d_IntersectionPoint(P, Uc1, Uc2, Trans1, Trans2,
    /// ReversedFlag). If ReversedFlag, the parameter/transition order swaps.
    pub fn new(
        p: DVec2,
        uc1: f64,
        uc2: f64,
        trans1: Transition,
        trans2: Transition,
        reversed_flag: bool,
    ) -> Self {
        if !reversed_flag {
            IntersectionPoint {
                pt: p,
                p1: uc1,
                p2: uc2,
                trans1,
                trans2,
            }
        } else {
            IntersectionPoint {
                pt: p,
                p1: uc2,
                p2: uc1,
                trans1: trans2,
                trans2: trans1,
            }
        }
    }

    pub fn set_values(
        &mut self,
        p: DVec2,
        uc1: f64,
        uc2: f64,
        trans1: Transition,
        trans2: Transition,
        reversed_flag: bool,
    ) {
        self.pt = p;
        if !reversed_flag {
            self.trans1 = trans1;
            self.trans2 = trans2;
            self.p1 = uc1;
            self.p2 = uc2;
        } else {
            self.trans1 = trans2;
            self.trans2 = trans1;
            self.p1 = uc2;
            self.p2 = uc1;
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

    pub fn transition_of_first(&self) -> &Transition {
        &self.trans1
    }

    pub fn transition_of_second(&self) -> &Transition {
        &self.trans2
    }
}

/// OCCT IntRes2d_IntersectionSegment — an intersection segment between two
/// 2D curves.
#[derive(Debug, Clone)]
pub struct IntersectionSegment {
    oppos: bool,
    first: bool,
    last: bool,
    ptfirst: IntersectionPoint,
    ptlast: IntersectionPoint,
}

impl IntersectionSegment {
    pub fn empty() -> Self {
        IntersectionSegment {
            oppos: false,
            first: false,
            last: false,
            ptfirst: IntersectionPoint::empty(),
            ptlast: IntersectionPoint::empty(),
        }
    }

    pub fn with_points(
        p1: &IntersectionPoint,
        p2: &IntersectionPoint,
        oppos: bool,
        reverse_flag: bool,
    ) -> Self {
        if !reverse_flag {
            IntersectionSegment {
                oppos,
                first: true,
                last: true,
                ptfirst: p1.clone(),
                ptlast: p2.clone(),
            }
        } else {
            IntersectionSegment {
                oppos,
                first: true,
                last: true,
                ptfirst: p2.clone(),
                ptlast: p1.clone(),
            }
        }
    }

    /// OCCT IntRes2d_IntersectionSegment(P, First, Oppos, ReverseFlag).
    pub fn with_one_point(
        p: &IntersectionPoint,
        first: bool,
        oppos: bool,
        reverse_flag: bool,
    ) -> Self {
        let (has_first, has_last, pt_first, pt_last) = if reverse_flag {
            if first {
                (false, true, IntersectionPoint::empty(), p.clone())
            } else {
                (true, false, p.clone(), IntersectionPoint::empty())
            }
        } else if first {
            (true, false, p.clone(), IntersectionPoint::empty())
        } else {
            (false, true, IntersectionPoint::empty(), p.clone())
        };
        IntersectionSegment {
            oppos,
            first: has_first,
            last: has_last,
            ptfirst: pt_first,
            ptlast: pt_last,
        }
    }

    /// Creates an infinite segment of intersection.
    pub fn infinite(oppos: bool) -> Self {
        IntersectionSegment {
            oppos,
            first: false,
            last: false,
            ptfirst: IntersectionPoint::empty(),
            ptlast: IntersectionPoint::empty(),
        }
    }

    pub fn is_opposite(&self) -> bool {
        self.oppos
    }

    pub fn has_first_point(&self) -> bool {
        self.first
    }

    pub fn first_point(&self) -> &IntersectionPoint {
        &self.ptfirst
    }

    pub fn has_last_point(&self) -> bool {
        self.last
    }

    pub fn last_point(&self) -> &IntersectionPoint {
        &self.ptlast
    }
}
