//! OCCT Law_S (TKGeomAlgo/Law) — 1:1 port of Law_S.hxx and Law_S.cxx
//! (whole file L25-52).
//!
//! Architecture mapping: `class Law_S : public Law_BSpFunc` is expressed by
//! composition on [`LawBSpFunc`]; every LawFunction method delegates to the
//! base subobject.

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::GeomAbsShape;

use super::law_bspline::LawBSpline;
use super::law_bsp_func::LawBSpFunc;
use super::law_function::{LawFunction, LawFunctionHandle};

/// OCCT Law_S — describes an "S" evolution law
/// (Law_S.hxx: `class Law_S : public Law_BSpFunc`).
#[derive(Debug, Clone)]
pub struct LawS {
    /// OCCT Law_BSpFunc base class subobject.
    pub base: LawBSpFunc,
}

impl LawS {
    /// OCCT Law_S::Law_S() (Law_S.cxx L25) — default (empty law).
    pub fn new() -> Self {
        LawS {
            base: LawBSpFunc::new(),
        }
    }

    /// OCCT Law_S::Set(Pdeb, Valdeb, Pfin, Valfin) (Law_S.cxx L27-30) —
    /// delegates to the 6-argument overload with null derivatives.
    pub fn set(&mut self, pdeb: f64, valdeb: f64, pfin: f64, valfin: f64) {
        self.set_with_tangents(pdeb, valdeb, 0.0, pfin, valfin, 0.0);
    }

    /// OCCT Law_S::Set(Pdeb, Valdeb, Ddeb, Pfin, Valfin, Dfin)
    /// (Law_S.cxx L32-52).
    pub fn set_with_tangents(
        &mut self,
        pdeb: f64,
        valdeb: f64,
        ddeb: f64,
        pfin: f64,
        valfin: f64,
        dfin: f64,
    ) {
        // OCCT 1-based arrays poles(1,4) / knots(1,2) / mults(1,2) — stored
        // 0-based in Rust.
        let mut poles = vec![0.0f64; 4];
        let knots = [pdeb, pfin];
        let mults = [4i32, 4];
        poles[0] = valdeb;
        poles[3] = valfin;
        let coe = (pfin - pdeb) / 3.0;
        poles[1] = valdeb + coe * ddeb;
        poles[2] = valfin - coe * dfin;

        let bs = Rc::new(RefCell::new(LawBSpline::new(&poles, &knots, &mults, 3, false)));
        // OCCT L51: SetCurve(new Law_BSpline(poles, knots, mults, 3));
        set_curve(&mut self.base, bs);
    }
}

/// OCCT Law_BSpFunc::SetCurve (Law_BSpFunc.cxx L396-401):
/// `curv = C; first = C->FirstParameter(); last = C->LastParameter();`
/// Written here as a helper because the rcad LawBSpFunc keeps its fields
/// private (law_bsp_func.rs predates this port); the handle-based OCCT
/// assignment maps to rebuilding the base value.
pub(crate) fn set_curve(base: &mut LawBSpFunc, c: Rc<RefCell<LawBSpline>>) {
    let (first, last) = {
        let curv = c.borrow();
        (curv.first_parameter(), curv.last_parameter())
    };
    *base = LawBSpFunc::with_curve(c, first, last);
}

// OCCT inheritance: every Law_Function virtual of Law_BSpFunc applies to
// Law_S; the impls below are pure delegation to the base subobject.
impl LawFunction for LawS {
    fn continuity(&self) -> GeomAbsShape {
        self.base.continuity()
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.base.nb_intervals(s)
    }

    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        self.base.intervals(t, s)
    }

    fn value(&mut self, x: f64) -> f64 {
        self.base.value(x)
    }

    fn d1(&mut self, x: f64, f: &mut f64, d: &mut f64) {
        self.base.d1(x, f, d)
    }

    fn d2(&mut self, x: f64, f: &mut f64, d: &mut f64, d2: &mut f64) {
        self.base.d2(x, f, d, d2)
    }

    fn trim(&self, pfirst: f64, plast: f64, tol: f64) -> LawFunctionHandle {
        self.base.trim(pfirst, plast, tol)
    }

    fn bounds(&self, pfirst: &mut f64, plast: &mut f64) {
        self.base.bounds(pfirst, plast)
    }
}
