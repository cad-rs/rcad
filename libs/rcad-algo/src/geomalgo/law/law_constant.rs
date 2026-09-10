//! OCCT Law_Constant (TKGeomAlgo/Law) — 1:1 port of Law_Constant.hxx
//! (L30-60) and Law_Constant.cxx (whole file L26-105).

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::GeomAbsShape;

use super::law_function::LawFunction;

/// OCCT Law_Constant — "Loi constante".
#[derive(Debug, Clone)]
pub struct LawConstant {
    radius: f64,
    first: f64,
    last: f64,
}

impl LawConstant {
    /// OCCT Law_Constant() (Law_Constant.cxx L31-36).
    pub fn new() -> Self {
        LawConstant {
            radius: 0.0,
            first: 0.0,
            last: 0.0,
        }
    }

    /// OCCT Set (L38-43) — set the radius and the range of the constant law.
    pub fn set(&mut self, radius: f64, pfirst: f64, plast: f64) {
        self.radius = radius;
        self.first = pfirst;
        self.last = plast;
    }
}

impl LawFunction for LawConstant {
    /// OCCT Continuity (L45-48) — returns GeomAbs_CN.
    fn continuity(&self) -> GeomAbsShape {
        GeomAbsShape::CN
    }

    /// OCCT NbIntervals (L50-53) — returns 1.
    fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        1
    }

    /// OCCT Intervals (L55-60).
    fn intervals(&self, t: &mut Vec<f64>, _s: GeomAbsShape) {
        let upper = t.len();
        t[0] = self.first;
        t[upper - 1] = self.last;
    }

    /// OCCT Value (L62-66).
    fn value(&mut self, _x: f64) -> f64 {
        self.radius
    }

    /// OCCT D1 (L68-72).
    fn d1(&mut self, _x: f64, f: &mut f64, d: &mut f64) {
        *f = self.radius;
        *d = 0.0;
    }

    /// OCCT D2 (L74-78).
    fn d2(&mut self, _x: f64, f: &mut f64, d: &mut f64, d2: &mut f64) {
        *f = self.radius;
        *d = 0.0;
        *d2 = 0.0;
    }

    /// OCCT Trim (L80-91).
    fn trim(&self, pfirst: f64, plast: f64, _tol: f64) -> super::law_function::LawFunctionHandle {
        let mut l = LawConstant::new();
        l.set(self.radius, pfirst, plast);
        Rc::new(RefCell::new(l))
    }

    /// OCCT Bounds (L93-98).
    fn bounds(&self, pfirst: &mut f64, plast: &mut f64) {
        *pfirst = self.first;
        *plast = self.last;
    }
}
