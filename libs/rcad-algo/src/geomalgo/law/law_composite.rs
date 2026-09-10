//! OCCT Law_Composite (TKGeomAlgo/Law) — 1:1 port of Law_Composite.hxx
//! (L30-70) and Law_Composite.cxx (whole file L26-247).
//!
//! Architecture mapping: `NCollection_List<handle<Law_Function>>` ->
//! `Vec<LawFunctionHandle>` (append-only use); the `curfunc` shared handle
//! member keeps its OCCT type (`Option<LawFunctionHandle>`).

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::el::in_period;
use rcad_kernel::math::GeomAbsShape;

use super::law_function::{LawFunction, LawFunctionHandle};

/// OCCT Law_Composite.
#[derive(Clone)]
pub struct LawComposite {
    first: f64,
    last: f64,
    curfunc: Option<LawFunctionHandle>,
    funclist: Vec<LawFunctionHandle>,
    periodic: bool,
    tfirst: f64,
    tlast: f64,
    ptol: f64,
}

impl LawComposite {
    /// OCCT Law_Composite() (Law_Composite.cxx L31-41).
    pub fn new() -> Self {
        LawComposite {
            first: -1.0e100,
            last: 1.0e100,
            curfunc: None,
            funclist: Vec::new(),
            periodic: false,
            tfirst: -1.0e100,
            tlast: 1.0e100,
            ptol: 0.0,
        }
    }

    /// OCCT Law_Composite(First, Last, Tol) (L43-53).
    pub fn with_tolerance(first: f64, last: f64, tol: f64) -> Self {
        LawComposite {
            first: -1.0e100,
            last: 1.0e100,
            curfunc: None,
            funclist: Vec::new(),
            periodic: false,
            tfirst: first,
            tlast: last,
            ptol: tol,
        }
    }

    /// OCCT ChangeElementaryLaw(W) (L202-207) — the elementary function of
    /// the composite used to compute at parameter W (shared handle).
    pub fn change_elementary_law(&mut self, w: f64) -> Option<&LawFunctionHandle> {
        let mut ww = w;
        self.prepare(&mut ww);
        self.curfunc.as_ref()
    }

    /// OCCT ChangeLaws() — the list of elementary laws.
    pub fn change_laws(&mut self) -> &mut Vec<LawFunctionHandle> {
        &mut self.funclist
    }

    /// OCCT IsPeriodic() (L209-212).
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// OCCT SetPeriodic() (L214-217).
    pub fn set_periodic(&mut self) {
        self.periodic = true;
    }

    /// OCCT Prepare (L152-200) — sets the current function.
    fn prepare(&mut self, w: &mut f64) {
        let mut f = 0.0;
        let mut l = 0.0;
        let eps;
        if *w - self.tfirst < self.tlast - *w {
            eps = self.ptol;
        } else {
            eps = -self.ptol;
        }
        if self.curfunc.is_none() {
            // OCCT: curfunc = funclist.Last(); curfunc->Bounds(f, last);
            //       curfunc = funclist.First(); curfunc->Bounds(first, l);
            self.curfunc = Some(self.funclist.last().unwrap().clone());
            let mut last_law = self.curfunc.as_ref().unwrap().borrow_mut();
            let mut f_tmp = 0.0;
            last_law.bounds(&mut f_tmp, &mut self.last);
            drop(last_law);
            self.curfunc = Some(self.funclist.first().unwrap().clone());
            let mut first_law = self.curfunc.as_ref().unwrap().borrow_mut();
            let mut l_tmp = 0.0;
            first_law.bounds(&mut self.first, &mut l_tmp);
            drop(first_law);
        }
        let mut wtest = *w + eps; // Decalage pour discriminer les noeuds
        if self.periodic {
            wtest = in_period(wtest, self.first, self.last);
            *w = wtest - eps;
        }
        {
            let cur = self.curfunc.as_ref().unwrap().borrow_mut();
            cur.bounds(&mut f, &mut l);
        }
        if f <= wtest && wtest <= l {
            return;
        }
        if *w <= self.first {
            self.curfunc = Some(self.funclist.first().unwrap().clone());
        } else if *w >= self.last {
            self.curfunc = Some(self.funclist.last().unwrap().clone());
        } else {
            for law in &self.funclist {
                self.curfunc = Some(law.clone());
                let cur = self.curfunc.as_ref().unwrap().borrow_mut();
                cur.bounds(&mut f, &mut l);
                drop(cur);
                if f <= wtest && wtest <= l {
                    return;
                }
            }
        }
    }
}

impl LawFunction for LawComposite {
    /// OCCT Continuity (L55-58).
    fn continuity(&self) -> GeomAbsShape {
        panic!("Law_Composite::Continuity()");
    }

    /// OCCT NbIntervals (L60-74).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let mut nbr_interval = 0usize;
        for func in &self.funclist {
            nbr_interval += func.borrow().nb_intervals(s);
        }
        nbr_interval
    }

    /// OCCT Intervals (L76-99).
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        let first = self.funclist.first().unwrap();
        let mut f_tmp = 0.0;
        let mut l_tmp = 0.0;
        first.borrow().bounds(&mut f_tmp, &mut l_tmp);
        t[0] = f_tmp;
        t[1] = l_tmp;
        let mut iglob = 2usize; // OCCT IGlob starts at 2 (1-based slot)
        for func in &self.funclist {
            let nb_index = func.borrow().nb_intervals(s) + 1;
            let mut loct: Vec<f64> = vec![0.0; nb_index];
            func.borrow().intervals(&mut loct, s);
            for iloc in 1..nb_index {
                // OCCT: for (Iloc = 2; Iloc <= nb_index; Iloc++, IGlob++)
                t[iglob] = loct[iloc];
                iglob += 1;
            }
        }
    }

    /// OCCT Value (L101-107).
    fn value(&mut self, x: f64) -> f64 {
        let mut w = x;
        self.prepare(&mut w);
        self.curfunc.as_ref().unwrap().borrow_mut().value(w)
    }

    /// OCCT D1 (L109-115).
    fn d1(&mut self, x: f64, f: &mut f64, d: &mut f64) {
        let mut w = x;
        self.prepare(&mut w);
        self.curfunc.as_ref().unwrap().borrow_mut().d1(w, f, d);
    }

    /// OCCT D2 (L117-123).
    fn d2(&mut self, x: f64, f: &mut f64, d: &mut f64, d2: &mut f64) {
        let mut w = x;
        self.prepare(&mut w);
        self.curfunc
            .as_ref()
            .unwrap()
            .borrow_mut()
            .d2(w, f, d, d2);
    }

    /// OCCT Trim (L125-135).
    fn trim(&self, pfirst: f64, plast: f64, tol: f64) -> LawFunctionHandle {
        let mut l = LawComposite::with_tolerance(pfirst, plast, tol);
        // l->ChangeLaws() = funclist — the handles are shared.
        l.change_laws().extend(self.funclist.iter().cloned());
        Rc::new(RefCell::new(l))
    }

    /// OCCT Bounds (L137-142).
    fn bounds(&self, pfirst: &mut f64, plast: &mut f64) {
        *pfirst = self.first;
        *plast = self.last;
    }
}
