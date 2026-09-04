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

    /// OCCT IntRes2d_Transition::SetValue(Tangent, Pos, Type) — sets an IN or
    /// OUT transition (IntRes2d_Transition.lxx L58-66).
    pub fn set_value_in_out(&mut self, tangent: bool, pos: Position, typetra: TypeTrans) {
        self.tangent = tangent;
        self.posit = pos;
        self.typetra = typetra;
    }

    /// OCCT IntRes2d_Transition::SetValue(Tangent, Pos, Situ, Oppos) — sets a
    /// TOUCH transition (IntRes2d_Transition.lxx L68-79).
    pub fn set_value_touch(&mut self, tangent: bool, pos: Position, situat: Situation, oppos: bool) {
        self.tangent = tangent;
        self.posit = pos;
        self.typetra = TypeTrans::Touch;
        self.situat = situat;
        self.oppos = oppos;
    }

    /// OCCT IntRes2d_Transition::SetValue(Pos) — sets an UNDECIDED transition
    /// (IntRes2d_Transition.lxx L81-86).
    pub fn set_value_undecided(&mut self, pos: Position) {
        self.posit = pos;
        self.typetra = TypeTrans::Undecided;
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

// OCCT IntRes2d_Intersection.cxx L24: PARAMEQUAL
const PARAMEQUAL_TOL: f64 = 1e-8;

fn paramequal(a: f64, b: f64) -> bool {
    (a - b).abs() < PARAMEQUAL_TOL
}

/// OCCT IntRes2d_Transition::TransitionEqual (IntRes2d_Intersection.cxx
/// L38-65) — equality of two transitions.
fn transition_equal(t1: &Transition, t2: &Transition) -> bool {
    if t1.position_on_curve() == t2.position_on_curve() {
        if t1.transition_type() == t2.transition_type() {
            if t1.transition_type() == TypeTrans::Touch {
                if t1.is_tangent() == t2.is_tangent() {
                    if t1.situation() == t2.situation() {
                        if t1.is_opposite() == t2.is_opposite() {
                            return true;
                        }
                    }
                }
            } else {
                return true;
            }
        }
    }
    false
}

/// OCCT IntRes2d_Intersection — root class of all the intersections between
/// two 2D curves (IntRes2d_Intersection.hxx/.lxx/.cxx). Holds the resulting
/// intersection points and segments plus the done/reverse flags.
#[derive(Debug, Clone)]
pub struct IntersectionBase {
    pub(crate) lpnt: Vec<IntersectionPoint>,
    pub(crate) lseg: Vec<IntersectionSegment>,
    pub(crate) done: bool,
    pub(crate) reverse: bool,
}

impl IntersectionBase {
    /// OCCT IntRes2d_Intersection() — empty constructor (done = reverse = false).
    pub fn new() -> Self {
        IntersectionBase {
            lpnt: Vec::new(),
            lseg: Vec::new(),
            done: false,
            reverse: false,
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsEmpty().
    pub fn is_empty(&self) -> bool {
        assert!(self.done, "StdFail_NotDone");
        self.lpnt.is_empty() && self.lseg.is_empty()
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        assert!(self.done, "StdFail_NotDone");
        self.lpnt.len()
    }

    /// OCCT Point(N) — 1-based index.
    pub fn point(&self, n: usize) -> &IntersectionPoint {
        assert!(self.done, "StdFail_NotDone");
        &self.lpnt[n - 1]
    }

    /// OCCT NbSegments().
    pub fn nb_segments(&self) -> usize {
        assert!(self.done, "StdFail_NotDone");
        self.lseg.len()
    }

    /// OCCT Segment(N) — 1-based index.
    pub fn segment(&self, n: usize) -> &IntersectionSegment {
        assert!(self.done, "StdFail_NotDone");
        &self.lseg[n - 1]
    }

    /// OCCT SetReversedParameters(flag).
    pub fn set_reversed_parameters(&mut self, flag: bool) {
        self.reverse = flag;
    }

    /// OCCT ReversedParameters().
    pub fn reversed_parameters(&self) -> bool {
        self.reverse
    }

    /// OCCT ResetFields() — clears the results and resets done.
    pub fn reset_fields(&mut self) {
        if self.done {
            self.lseg.clear();
            self.lpnt.clear();
            self.done = false;
        }
    }

    /// OCCT Insert(Pnt) (IntRes2d_Intersection.cxx L67-113) — inserts the point
    /// sorted by ParamOnFirst, skipping an exact duplicate (parameter + both
    /// transitions equal).
    pub fn insert(&mut self, pnt: &IntersectionPoint) {
        let n = self.lpnt.len();
        if n == 0 {
            self.lpnt.push(pnt.clone());
            return;
        }
        let u = pnt.param_on_first();
        let mut i = 0usize;
        let mut b = n + 1;
        while i < n {
            let ui = self.lpnt[i].param_on_first();
            if ui >= u {
                b = i + 1;
                i = n;
            }
            if paramequal(ui, u) {
                if paramequal(pnt.param_on_second(), self.lpnt[i].param_on_second()) {
                    if transition_equal(pnt.transition_of_first(), self.lpnt[i].transition_of_first())
                        && transition_equal(
                            pnt.transition_of_second(),
                            self.lpnt[i].transition_of_second(),
                        )
                    {
                        b = 0;
                        i = n;
                    }
                }
            }
            i += 1;
        }
        if b > n {
            self.lpnt.push(pnt.clone());
        } else if b > 0 {
            self.lpnt.insert(b - 1, pnt.clone());
        }
    }

    /// OCCT Append(Pnt) — appends a point.
    pub fn append_point(&mut self, pnt: &IntersectionPoint) {
        self.lpnt.push(pnt.clone());
    }

    /// OCCT Append(Seg) — appends a segment.
    pub fn append_segment(&mut self, seg: &IntersectionSegment) {
        self.lseg.push(seg.clone());
    }

    /// OCCT SetValues(Other) (IntRes2d_Intersection.cxx L115-142) — copies the
    /// results of another intersection.
    pub fn set_values(&mut self, other: &IntersectionBase) {
        if other.done {
            self.lseg.clear();
            self.lpnt.clear();
            for p in other.lpnt.iter() {
                self.lpnt.push(p.clone());
            }
            for s in other.lseg.iter() {
                self.lseg.push(s.clone());
            }
            self.done = true;
        } else {
            self.done = false;
        }
    }

    /// OCCT Append(Other, FirstParam1, LastParam1, FirstParam2, LastParam2)
    /// (IntRes2d_Intersection.cxx L182-390) — merges the results of another
    /// intersection restricted to the parameter windows, joining
    /// collinear-continuation segments and dropping points interior to a
    /// merged segment.
    pub fn append_intersector(
        &mut self,
        other: &IntersectionBase,
        first_param1: f64,
        last_param1: f64,
        first_param2: f64,
        last_param2: f64,
    ) {
        if other.done {
            // -- Verification of the Position of the IntersectionPoints.
            let n = other.lpnt.len();
            for i in 1..=n {
                let p = &other.lpnt[i - 1];
                let p_param_on_first = p.param_on_first();
                let p_param_on_second = p.param_on_second();
                let mut t1 = p.transition_of_first().clone();
                let mut t2 = p.transition_of_second().clone();
                let pt = p.value();

                internal_verify_position(
                    &mut t1,
                    &mut t2,
                    p_param_on_first,
                    p_param_on_second,
                    first_param1,
                    last_param1,
                    first_param2,
                    last_param2,
                );

                self.insert(&IntersectionPoint::new(
                    pt,
                    p_param_on_first,
                    p_param_on_second,
                    t1,
                    t2,
                    false,
                ));
            }

            //--------------------------------------------------
            //-- IntersectionSegment
            //-- (we assume that a composite curve is always bounded)
            //-- (a segment has always a FirstPoint and a LastPoint)
            //--------------------------------------------------
            let n = other.lseg.len();
            let (mut seg_modif_p1first, mut seg_modif_p1second) = (0.0, 0.0);
            let (mut seg_modif_p2first, mut seg_modif_p2second) = (0.0, 0.0);

            for i in 1..=n {
                let p1 = other.lseg[i - 1].first_point().clone();
                let p1_param_on_first = p1.param_on_first();
                let p1_param_on_second = p1.param_on_second();
                let mut p1_t1 = p1.transition_of_first().clone();
                let mut p1_t2 = p1.transition_of_second().clone();
                let p1_pt = p1.value();

                internal_verify_position(
                    &mut p1_t1,
                    &mut p1_t2,
                    p1_param_on_first,
                    p1_param_on_second,
                    first_param1,
                    last_param1,
                    first_param2,
                    last_param2,
                );

                let p2 = other.lseg[i - 1].last_point().clone();
                let p2_param_on_first = p2.param_on_first();
                let p2_param_on_second = p2.param_on_second();
                let mut p2_t1 = p2.transition_of_first().clone();
                let mut p2_t2 = p2.transition_of_second().clone();
                let p2_pt = p2.value();

                let oppos = other.lseg[i - 1].is_opposite();

                internal_verify_position(
                    &mut p2_t1,
                    &mut p2_t2,
                    p2_param_on_first,
                    p2_param_on_second,
                    first_param1,
                    last_param1,
                    first_param2,
                    last_param2,
                );

                //-- Loop on the previous segments.
                let an = self.lseg.len();
                let mut not_yet_modified = true;
                let mut j = 1usize;
                while (j <= an) && not_yet_modified {
                    let anp1 = self.lseg[j - 1].first_point().clone();
                    let anp1_param_on_first = anp1.param_on_first();
                    let anp1_param_on_second = anp1.param_on_second();
                    let anp2 = self.lseg[j - 1].last_point().clone();
                    let anp2_param_on_first = anp2.param_on_first();
                    let anp2_param_on_second = anp2.param_on_second();

                    if oppos == self.lseg[j - 1].is_opposite() {
                        //---------------------------------------------------------------
                        //--    AnP1---------AnP2
                        //--                  P1-------------P2
                        //--
                        if paramequal(p1_param_on_first, anp2_param_on_first)
                            && paramequal(p1_param_on_second, anp2_param_on_second)
                        {
                            not_yet_modified = false;
                            self.lseg[j - 1] =
                                IntersectionSegment::with_points(&anp1, &p2, oppos, false);
                            seg_modif_p1first = anp1_param_on_first;
                            seg_modif_p1second = anp1_param_on_second;
                            seg_modif_p2first = p2_param_on_first;
                            seg_modif_p2second = p2_param_on_second;
                        }
                        //---------------------------------------------------------------
                        //--                                AnP1---------AnP2
                        //--                  P1-------------P2
                        //--
                        else if paramequal(p2_param_on_first, anp1_param_on_first)
                            && paramequal(p2_param_on_second, anp1_param_on_second)
                        {
                            not_yet_modified = false;
                            self.lseg[j - 1] =
                                IntersectionSegment::with_points(&p1, &anp2, oppos, false);
                            seg_modif_p1first = p1_param_on_first;
                            seg_modif_p1second = p1_param_on_second;
                            seg_modif_p2first = anp2_param_on_first;
                            seg_modif_p2second = anp2_param_on_second;
                        }
                        //---------------------------------------------------------------
                        //--    AnP2---------AnP1
                        //--                  P1-------------P2
                        //--
                        if paramequal(p1_param_on_first, anp1_param_on_first)
                            && paramequal(p1_param_on_second, anp1_param_on_second)
                        {
                            not_yet_modified = false;
                            self.lseg[j - 1] =
                                IntersectionSegment::with_points(&anp2, &p2, oppos, false);
                            seg_modif_p1first = p2_param_on_first;
                            seg_modif_p1second = p2_param_on_second;
                            seg_modif_p2first = anp2_param_on_first;
                            seg_modif_p2second = anp2_param_on_second;
                        }
                        //---------------------------------------------------------------
                        //--                                AnP2---------AnP1
                        //--                  P1-------------P2
                        //--
                        else if paramequal(p2_param_on_first, anp2_param_on_first)
                            && paramequal(p2_param_on_second, anp2_param_on_second)
                        {
                            not_yet_modified = false;
                            self.lseg[j - 1] =
                                IntersectionSegment::with_points(&p1, &anp1, oppos, false);
                            seg_modif_p1first = p1_param_on_first;
                            seg_modif_p1second = p1_param_on_second;
                            seg_modif_p2first = anp1_param_on_first;
                            seg_modif_p2second = anp1_param_on_second;
                        }
                    }
                    j += 1;
                }
                if not_yet_modified {
                    self.append_segment(&IntersectionSegment::with_points(
                        &IntersectionPoint::new(
                            p1_pt,
                            p1_param_on_first,
                            p1_param_on_second,
                            p1_t1,
                            p1_t2,
                            false,
                        ),
                        &IntersectionPoint::new(
                            p2_pt,
                            p2_param_on_first,
                            p2_param_on_second,
                            p2_t1,
                            p2_t2,
                            false,
                        ),
                        oppos,
                        false,
                    ));
                } else {
                    //--------------------------------------------------------------
                    //-- Are some Existing Points in this segment ?
                    //--------------------------------------------------------------
                    let mut rnbpts = self.lpnt.len() as i64;
                    let mut rp: i64 = 1;
                    while (rp <= rnbpts) && (rp >= 1) {
                        let pon_first = self.lpnt[(rp - 1) as usize].param_on_first();
                        let pon_second = self.lpnt[(rp - 1) as usize].param_on_second();

                        if ((pon_first >= seg_modif_p1first && pon_first <= seg_modif_p2first)
                            || (pon_first <= seg_modif_p1first && pon_first >= seg_modif_p2first))
                            && ((pon_second >= seg_modif_p1second
                                && pon_second <= seg_modif_p2second)
                                || (pon_second <= seg_modif_p1second
                                    && pon_second >= seg_modif_p2second))
                        {
                            self.lpnt.remove((rp - 1) as usize);
                            rp -= 1;
                            rnbpts -= 1;
                        }
                        rp += 1;
                    }
                }
            }
            //--------------------------------------------------
            //-- Remove some Points ?
            //-- Example : Points which lie in a segment.
            //--------------------------------------------------

            self.done = true;
        } else {
            self.done = false;
        }
    }
}

/// OCCT static InternalVerifyPosition (IntRes2d_Intersection.cxx L431-465).
#[allow(clippy::too_many_arguments)]
fn internal_verify_position(
    t1: &mut Transition,
    t2: &mut Transition,
    p_param_on_first: f64,
    p_param_on_second: f64,
    first_param1: f64,
    last_param1: f64,
    first_param2: f64,
    last_param2: f64,
) {
    if t1.position_on_curve() != Position::Middle
        && !(paramequal(p_param_on_first, first_param1) || paramequal(p_param_on_first, last_param1))
    {
        if (p_param_on_first > first_param1) && (p_param_on_first < last_param1) {
            t1.set_position(Position::Middle);
        }
    }
    if t2.position_on_curve() != Position::Middle
        && !(paramequal(p_param_on_second, first_param2)
            || paramequal(p_param_on_second, last_param2))
    {
        if (p_param_on_second > first_param2) && (p_param_on_second < last_param2) {
            t2.set_position(Position::Middle);
        }
    }
}

impl Default for IntersectionBase {
    fn default() -> Self {
        IntersectionBase::new()
    }
}
