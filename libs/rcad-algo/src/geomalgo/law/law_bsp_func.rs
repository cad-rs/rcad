//! OCCT Law_BSpFunc (TKGeomAlgo/Law) — 1:1 port of Law_BSpFunc.cxx (whole
//! file L26-401) + .hxx members.  The `PosTol` constant is
//! Precision::PConfusion() / 2.

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::bspl_lib::locate_parameter_knots_mults;

use super::law_bspline::LawBSpline;
use super::law_bspline_knot_splitting::LawBSplineKnotSplitting;
use super::law_function::{LawFunction, LawFunctionHandle};

// OCCT: #define PosTol Precision::PConfusion() / 2
const POS_TOL: f64 = 1e-12 / 2.0;

/// OCCT Law_BSpFunc — a Law_Function built on a Law_BSpline.
#[derive(Debug, Clone)]
pub struct LawBSpFunc {
    curv: Option<Rc<RefCell<LawBSpline>>>,
    first: f64,
    last: f64,
}

impl LawBSpFunc {
    /// OCCT Law_BSpFunc() (L36-40).
    pub fn new() -> Self {
        LawBSpFunc {
            curv: None,
            first: 0.0,
            last: 0.0,
        }
    }

    /// OCCT Law_BSpFunc(C, First, Last) (L42-48).
    pub fn with_curve(curv: Rc<RefCell<LawBSpline>>, first: f64, last: f64) -> Self {
        LawBSpFunc {
            curv: Some(curv),
            first,
            last,
        }
    }

    /// OCCT Curve() — the underlying Law_BSpline handle.
    pub fn curve(&self) -> Option<&Rc<RefCell<LawBSpline>>> {
        self.curv.as_ref()
    }
}

impl LawFunction for LawBSpFunc {
    /// OCCT Continuity (L50-53) — `return curv->Continuity()`.  The rcad
    /// LawBSpline continuity is a SmoothShape (C0..C3/CN; the G1/G2 values
    /// are never produced by updateKnots), mapped back to GeomAbsShape.
    fn continuity(&self) -> GeomAbsShape {
        use super::law_bspline::SmoothShape;
        let curv = self.curv.as_ref().unwrap().borrow();
        match curv.continuity() {
            SmoothShape::C0 => GeomAbsShape::C0,
            SmoothShape::C1 => GeomAbsShape::C1,
            SmoothShape::C2 => GeomAbsShape::C2,
            SmoothShape::C3 => GeomAbsShape::C3,
            SmoothShape::CN => GeomAbsShape::CN,
            _ => unreachable!(),
        }
    }

    /// OCCT NbIntervals (L55-136).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let mut my_nb_intervals = 1usize;
        let curv = self.curv.as_ref().unwrap().borrow();
        let cont_of_curve = curve_continuity_rank(&curv);
        let s_rank = shape_rank(s);
        if s_rank > cont_of_curve {
            match s {
                GeomAbsShape::C0 => {
                    my_nb_intervals = 1;
                }
                GeomAbsShape::C1
                | GeomAbsShape::C2
                | GeomAbsShape::C3
                | GeomAbsShape::CN => {
                    let cont = match s {
                        GeomAbsShape::C1 => 1,
                        GeomAbsShape::C2 => 2,
                        GeomAbsShape::C3 => 3,
                        _ => curv.degree(),
                    };
                    let convector = LawBSplineKnotSplitting::new(&curv, cont as i32);
                    let nb_int = convector.nb_splits() - 1;
                    let mut inter = vec![0i32; nb_int + 1];
                    convector.splitting(&mut inter);
                    let nb = curv.nb_knots();
                    let mut index1 = 0i32;
                    let mut index2 = 0i32;
                    let mut new_first = 0.0;
                    let mut new_last = 0.0;
                    let tk = curv.knots().clone();
                    let tm = curv.multiplicities().clone();
                    locate_parameter_knots_mults(
                        curv.degree(),
                        &tk,
                        &tm,
                        self.first,
                        curv.is_periodic(),
                        1,
                        nb as i32,
                        &mut index1,
                        &mut new_first,
                    );
                    locate_parameter_knots_mults(
                        curv.degree(),
                        &tk,
                        &tm,
                        self.last,
                        curv.is_periodic(),
                        1,
                        nb as i32,
                        &mut index2,
                        &mut new_last,
                    );
                    let mut index1 = index1 as usize;
                    let mut index2 = index2 as usize;
                    if (new_first - tk[index1]).abs() < 1e-12 {
                        index1 += 1;
                    }
                    if new_last - tk[index2 - 1] > 1e-12 {
                        index2 += 1;
                    }
                    my_nb_intervals = 1;
                    for i in 1..=nb_int {
                        if (inter[i - 1] as usize) > index1 && (inter[i - 1] as usize) < index2 {
                            my_nb_intervals += 1;
                        }
                    }
                }
            }
        }
        my_nb_intervals
    }

    /// OCCT Intervals (L138-231).
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        let mut my_nb_intervals = 1usize;
        let curv = self.curv.as_ref().unwrap().borrow();
        let cont_of_curve = curve_continuity_rank(&curv);
        let s_rank = shape_rank(s);
        if s_rank > cont_of_curve {
            match s {
                GeomAbsShape::C0 => {
                    my_nb_intervals = 1;
                }
                GeomAbsShape::C1
                | GeomAbsShape::C2
                | GeomAbsShape::C3
                | GeomAbsShape::CN => {
                    let cont = match s {
                        GeomAbsShape::C1 => 1,
                        GeomAbsShape::C2 => 2,
                        GeomAbsShape::C3 => 3,
                        _ => curv.degree(),
                    };
                    let convector = LawBSplineKnotSplitting::new(&curv, cont as i32);
                    let nb_int = convector.nb_splits() - 1;
                    let mut inter = vec![0i32; nb_int + 1];
                    convector.splitting(&mut inter);
                    let nb = curv.nb_knots();
                    let mut index1 = 0i32;
                    let mut index2 = 0i32;
                    let mut new_first = 0.0;
                    let mut new_last = 0.0;
                    let tk = curv.knots().clone();
                    let tm = curv.multiplicities().clone();
                    locate_parameter_knots_mults(
                        curv.degree(),
                        &tk,
                        &tm,
                        self.first,
                        curv.is_periodic(),
                        1,
                        nb as i32,
                        &mut index1,
                        &mut new_first,
                    );
                    locate_parameter_knots_mults(
                        curv.degree(),
                        &tk,
                        &tm,
                        self.last,
                        curv.is_periodic(),
                        1,
                        nb as i32,
                        &mut index2,
                        &mut new_last,
                    );
                    let mut index1 = index1 as usize;
                    let mut index2 = index2 as usize;
                    if (new_first - tk[index1]).abs() < 1e-12 {
                        index1 += 1;
                    }
                    if new_last - tk[index2 - 1] > 1e-12 {
                        index2 += 1;
                    }
                    inter[0] = index1 as i32;
                    my_nb_intervals = 1;
                    for i in 1..=nb_int {
                        if (inter[i - 1] as usize) > index1 && (inter[i - 1] as usize) < index2 {
                            my_nb_intervals += 1;
                            inter[my_nb_intervals - 1] = inter[i - 1];
                        }
                    }
                    inter[my_nb_intervals] = index2 as i32;
                    for i in 1..=my_nb_intervals + 1 {
                        t[i - 1] = tk[(inter[i - 1] - 1) as usize];
                    }
                }
            }
        }
        t[0] = self.first;
        t[my_nb_intervals] = self.last;
    }

    /// OCCT Value (L233-269).
    fn value(&mut self, x: f64) -> f64 {
        if x == self.first || x == self.last {
            let (mut ideb, mut ifin) = (0i32, 0i32);
            let curv = self.curv.as_ref().unwrap().borrow();
            if x == self.first {
                curv.locate_u(self.first, POS_TOL, &mut ideb, &mut ifin, true);
                if ideb < 1 {
                    ideb = 1;
                }
                if ideb >= ifin {
                    ifin = ideb + 1;
                }
            }
            if x == self.last {
                curv.locate_u(self.last, POS_TOL, &mut ideb, &mut ifin, true);
                if ifin > curv.nb_knots() as i32 {
                    ifin = curv.nb_knots() as i32;
                }
                if ideb >= ifin {
                    ideb = ifin - 1;
                }
            }
            curv.local_value(x, ideb, ifin)
        } else {
            self.curv.as_ref().unwrap().borrow().value(x)
        }
    }

    /// OCCT D1 (L271-303).
    fn d1(&mut self, x: f64, f: &mut f64, d: &mut f64) {
        if x == self.first || x == self.last {
            let (mut ideb, mut ifin) = (0i32, 0i32);
            let curv = self.curv.as_ref().unwrap().borrow();
            if x == self.first {
                curv.locate_u(self.first, POS_TOL, &mut ideb, &mut ifin, true);
                if ideb < 1 {
                    ideb = 1;
                }
                if ideb >= ifin {
                    ifin = ideb + 1;
                }
            }
            if x == self.last {
                curv.locate_u(self.last, POS_TOL, &mut ideb, &mut ifin, true);
                if ifin > curv.nb_knots() as i32 {
                    ifin = curv.nb_knots() as i32;
                }
                if ideb >= ifin {
                    ideb = ifin - 1;
                }
            }
            curv.local_d1(x, ideb, ifin, f, d);
        } else {
            self.curv.as_ref().unwrap().borrow_mut().d1(x, f, d);
        }
    }

    /// OCCT D2 (L305-337).
    fn d2(&mut self, x: f64, f: &mut f64, d: &mut f64, d2: &mut f64) {
        if x == self.first || x == self.last {
            let (mut ideb, mut ifin) = (0i32, 0i32);
            let curv = self.curv.as_ref().unwrap().borrow();
            if x == self.first {
                curv.locate_u(self.first, POS_TOL, &mut ideb, &mut ifin, true);
                if ideb < 1 {
                    ideb = 1;
                }
                if ideb >= ifin {
                    ifin = ideb + 1;
                }
            }
            if x == self.last {
                curv.locate_u(self.last, POS_TOL, &mut ideb, &mut ifin, true);
                if ifin > curv.nb_knots() as i32 {
                    ifin = curv.nb_knots() as i32;
                }
                if ideb >= ifin {
                    ideb = ifin - 1;
                }
            }
            curv.local_d2(x, ideb, ifin, f, d, d2);
        } else {
            self.curv.as_ref().unwrap().borrow_mut().d2(x, f, d, d2);
        }
    }

    /// OCCT Trim (L339-347).
    fn trim(&self, pfirst: f64, plast: f64, _tol: f64) -> LawFunctionHandle {
        let l = LawBSpFunc::with_curve(self.curv.as_ref().unwrap().clone(), pfirst, plast);
        Rc::new(RefCell::new(l))
    }

    /// OCCT Bounds (L349-353).
    fn bounds(&self, pfirst: &mut f64, plast: &mut f64) {
        *pfirst = self.first;
        *plast = self.last;
    }
}

// The LawBSpline continuity is a SmoothShape; rank it against GeomAbsShape
// the way the OCCT `S > Continuity()` comparison does (GeomAbs ordering:
// C0=0, G1=1, C1=2, G2=3, C2=4, C3=5, CN=6).
fn curve_continuity_rank(curv: &LawBSpline) -> i32 {
    use super::law_bspline::SmoothShape;
    match curv.continuity() {
        SmoothShape::C0 => 0,
        SmoothShape::G1 => 1,
        SmoothShape::C1 => 2,
        SmoothShape::G2 => 3,
        SmoothShape::C2 => 4,
        SmoothShape::C3 => 5,
        SmoothShape::CN => 6,
    }
}

fn shape_rank(s: GeomAbsShape) -> i32 {
    match s {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::C1 => 2,
        GeomAbsShape::C2 => 4,
        GeomAbsShape::C3 => 5,
        GeomAbsShape::CN => 6,
    }
}
