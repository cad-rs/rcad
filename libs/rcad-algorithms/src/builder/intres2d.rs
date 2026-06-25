use glam::DVec2;

#[derive(Debug, Clone)]
pub struct IntRes2dDomain {
    pub status: i32, pub first_param: f64, pub last_param: f64,
    pub first_tol: f64, pub last_tol: f64,
    pub first_point: DVec2, pub last_point: DVec2,
    pub period_first: f64, pub period_last: f64,
}
impl IntRes2dDomain {
    pub fn new() -> Self { IntRes2dDomain { status: 0, first_param: 0.0, last_param: 0.0, first_tol: 0.0, last_tol: 0.0, first_point: DVec2::ZERO, last_point: DVec2::ZERO, period_first: 0.0, period_last: 0.0 } }
    pub fn new_bounded(p1: DVec2, par1: f64, tol1: f64, p2: DVec2, par2: f64, tol2: f64) -> Self { let mut d = Self::new(); d.set_values_bounded(p1, par1, tol1, p2, par2, tol2); d }
    pub fn set_values_bounded(&mut self, p1: DVec2, par1: f64, tol1: f64, p2: DVec2, par2: f64, tol2: f64) { self.status = 3; self.period_first = 0.0; self.period_last = 0.0; self.first_param = par1; self.first_point = p1; self.first_tol = tol1; self.last_param = par2; self.last_point = p2; self.last_tol = tol2; }
    pub fn set_equivalent_parameters(&mut self, zero: f64, period: f64) { debug_assert!((self.status & 3) == 3); self.status |= 4; self.period_first = zero; self.period_last = period; }
    pub fn has_first_point(&self) -> bool { (self.status & 1) != 0 }
    pub fn first_parameter(&self) -> f64 { debug_assert!(self.has_first_point()); self.first_param }
    pub fn has_last_point(&self) -> bool { (self.status & 2) != 0 }
    pub fn last_parameter(&self) -> f64 { debug_assert!(self.has_last_point()); self.last_param }
    pub fn is_closed(&self) -> bool { (self.status & 4) != 0 }
    pub fn equivalent_parameters(&self) -> (f64, f64) { (self.period_first, self.period_last) }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntRes2dPosition { Head, Middle, End }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntRes2dTypeTrans { In, Out, Touch, Undecided }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntRes2dSituation { Inside, Outside, Unknown }

#[derive(Debug, Clone)]
pub struct IntRes2dTransition { pub tangent: bool, pub posit: IntRes2dPosition, pub typetra: IntRes2dTypeTrans, pub situat: IntRes2dSituation, pub oppos: bool }
impl IntRes2dTransition {
    pub fn new_in_out(tangent: bool, pos: IntRes2dPosition, typ: IntRes2dTypeTrans) -> Self { IntRes2dTransition { tangent, posit: pos, typetra: typ, situat: IntRes2dSituation::Unknown, oppos: false } }
    pub fn new_touch(tangent: bool, pos: IntRes2dPosition, situ: IntRes2dSituation, oppos: bool) -> Self { IntRes2dTransition { tangent, posit: pos, typetra: IntRes2dTypeTrans::Touch, situat: situ, oppos } }
    pub fn new_undecided(pos: IntRes2dPosition) -> Self { IntRes2dTransition { tangent: true, posit: pos, typetra: IntRes2dTypeTrans::Undecided, situat: IntRes2dSituation::Unknown, oppos: false } }
    pub fn transition_type(&self) -> IntRes2dTypeTrans { self.typetra }
    pub fn is_tangent(&self) -> bool { self.tangent }
    pub fn situation(&self) -> IntRes2dSituation { self.situat }
    pub fn is_opposite(&self) -> bool { self.oppos }
}

#[derive(Debug, Clone)]
pub struct IntRes2dIntersectionPoint { pub pt: DVec2, pub p1: f64, pub p2: f64, pub trans1: IntRes2dTransition, pub trans2: IntRes2dTransition }
impl IntRes2dIntersectionPoint {
    pub fn new(p: DVec2, u1: f64, u2: f64, t1: IntRes2dTransition, t2: IntRes2dTransition, reversed: bool) -> Self {
        if reversed { IntRes2dIntersectionPoint { pt: p, p1: u2, p2: u1, trans1: t2, trans2: t1 } }
        else { IntRes2dIntersectionPoint { pt: p, p1: u1, p2: u2, trans1: t1, trans2: t2 } }
    }
    pub fn value(&self) -> DVec2 { self.pt }
    pub fn param_on_first(&self) -> f64 { self.p1 }
    pub fn param_on_second(&self) -> f64 { self.p2 }
}
