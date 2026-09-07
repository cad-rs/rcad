//! OCCT Law_Linear (TKGeomAlgo/Law) — 1:1 port of Law_Linear.hxx
//! (L27-...) and Law_Linear.cxx (whole file L25-104).

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::GeomAbsShape;

use super::law_function::{LawFunction, LawFunctionHandle};

/// OCCT Law_Linear — describes a linear evolution law
/// (Law_Linear.hxx: `class Law_Linear : public Law_Function`).
#[derive(Debug, Clone)]
pub struct LawLinear {
    valdeb: f64,
    valfin: f64,
    pdeb: f64,
    pfin: f64,
}

impl LawLinear {
    /// OCCT Law_Linear::Law_Linear() (Law_Linear.cxx L25-31).
    pub fn new() -> Self {
        LawLinear {
            valdeb: 0.0,
            valfin: 0.0,
            pdeb: 0.0,
            pfin: 0.0,
        }
    }

    /// OCCT Law_Linear::Set (Law_Linear.cxx L33-39).
    pub fn set(&mut self, pdeb: f64, valdeb: f64, pfin: f64, valfin: f64) {
        self.pdeb = pdeb;
        self.pfin = pfin;
        self.valdeb = valdeb;
        self.valfin = valfin;
    }
}

impl LawFunction for LawLinear {
    /// OCCT Law_Linear::Continuity (Law_Linear.cxx L43-46) — returns GeomAbs_CN.
    fn continuity(&self) -> GeomAbsShape {
        GeomAbsShape::CN
    }

    /// OCCT Law_Linear::NbIntervals (Law_Linear.cxx L51-54) — returns 1.
    fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        1
    }

    /// OCCT Law_Linear::Intervals (Law_Linear.cxx L58-64).
    fn intervals(&self, t: &mut Vec<f64>, _s: GeomAbsShape) {
        let lower = 0usize; // OCCT: T.Lower()
        let upper = t.len() - 1; // OCCT: T.Upper()
        t[lower] = self.pdeb;
        t[upper] = self.pfin;
    }

    /// OCCT Law_Linear::Value (Law_Linear.cxx L66-69).
    fn value(&mut self, x: f64) -> f64 {
        ((x - self.pdeb) * self.valfin + (self.pfin - x) * self.valdeb) / (self.pfin - self.pdeb)
    }

    /// OCCT Law_Linear::D1 (Law_Linear.cxx L71-75).
    fn d1(&mut self, x: f64, f: &mut f64, d: &mut f64) {
        *f = ((x - self.pdeb) * self.valfin + (self.pfin - x) * self.valdeb)
            / (self.pfin - self.pdeb);
        *d = (self.valfin - self.valdeb) / (self.pfin - self.pdeb);
    }

    /// OCCT Law_Linear::D2 (Law_Linear.cxx L77-82).
    fn d2(&mut self, x: f64, f: &mut f64, d: &mut f64, d2: &mut f64) {
        *f = ((x - self.pdeb) * self.valfin + (self.pfin - x) * self.valdeb)
            / (self.pfin - self.pdeb);
        *d = (self.valfin - self.valdeb) / (self.pfin - self.pdeb);
        *d2 = 0.0;
    }

    /// OCCT Law_Linear::Trim (Law_Linear.cxx L86-98).
    fn trim(&self, pfirst: f64, plast: f64, _tol: f64) -> LawFunctionHandle {
        let mut l = LawLinear::new();
        let vdeb = ((pfirst - self.pdeb) * self.valfin + (self.pfin - pfirst) * self.valdeb)
            / (self.pfin - self.pdeb);
        let vfin = ((plast - self.pdeb) * self.valfin + (self.pfin - plast) * self.valdeb)
            / (self.pfin - self.pdeb);
        l.set(pfirst, vdeb, plast, vfin);
        Rc::new(RefCell::new(l))
    }

    /// OCCT Law_Linear::Bounds (Law_Linear.cxx L100-104).
    fn bounds(&self, pfirst: &mut f64, plast: &mut f64) {
        *pfirst = self.pdeb;
        *plast = self.pfin;
    }
}
