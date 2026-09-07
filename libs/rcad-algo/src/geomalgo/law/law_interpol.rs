//! OCCT Law_Interpol (TKGeomAlgo/Law) — 1:1 port of Law_Interpol.hxx and
//! Law_Interpol.cxx (whole file L32-181).
//!
//! Architecture mapping: `class Law_Interpol : public Law_BSpFunc` is
//! expressed by composition on [`LawBSpFunc`]; `NCollection_Array1<gp_Pnt2d>`
//! maps to a `&[DVec2]` slice (0-based stand-in for the OCCT 1-based array).

use glam::DVec2;

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::math::GeomAbsShape;

use super::law_bsp_func::LawBSpFunc;
use super::law_function::{LawFunction, LawFunctionHandle};
use super::law_interpolate::LawInterpolate;

/// OCCT Law_Interpol — provides an evolution law that interpolates a set
/// of parameter and value pairs (wi, radi)
/// (Law_Interpol.hxx: `class Law_Interpol : public Law_BSpFunc`).
#[derive(Debug, Clone)]
pub struct LawInterpol {
    /// OCCT Law_BSpFunc base class subobject.
    pub base: LawBSpFunc,
}

impl LawInterpol {
    /// OCCT Law_Interpol::Law_Interpol() (Law_Interpol.cxx L32) —
    /// constructs an empty interpolative evolution law.
    pub fn new() -> Self {
        LawInterpol {
            base: LawBSpFunc::new(),
        }
    }

    /// OCCT Law_Interpol::Set(ParAndRad, Periodic) (Law_Interpol.cxx L36-67).
    pub fn set(&mut self, par_and_rad: &[DVec2], periodic: bool) {
        let l = 0usize; // OCCT: ParAndRad.Lower()
        let nbp = par_and_rad.len(); // OCCT: ParAndRad.Length()

        // OCCT: par = new NCollection_HArray1<double>(1, nbp);
        let mut par = vec![0.0f64; nbp];
        // OCCT: rad has nbp-1 entries when Periodic, nbp otherwise.
        let rad_len = if periodic { nbp - 1 } else { nbp };
        let mut rad = vec![0.0f64; rad_len];
        for i in 1..=nbp {
            let coord = par_and_rad[l + i - 1];
            let (x, y) = (coord.x, coord.y);
            par[i - 1] = x;
            if !periodic || i != nbp {
                rad[i - 1] = y;
            }
        }
        let mut inter = LawInterpolate::with_parameters(rad, par, periodic, CONFUSION);
        inter.perform();
        let curve = inter.curve();
        super::law_s::set_curve(&mut self.base, curve);
    }

    /// OCCT Law_Interpol::SetInRelative(ParAndRad, Ud, Uf, Periodic)
    /// (Law_Interpol.cxx L71-105).
    pub fn set_in_relative(
        &mut self,
        par_and_rad: &[DVec2],
        ud: f64,
        uf: f64,
        periodic: bool,
    ) {
        let l = 0usize; // OCCT: ParAndRad.Lower()
        let u = par_and_rad.len() - 1; // OCCT: ParAndRad.Upper()
        let wd = par_and_rad[l].x;
        let wf = par_and_rad[u].x;
        let nbp = par_and_rad.len(); // OCCT: ParAndRad.Length()

        let mut par = vec![0.0f64; nbp];
        let rad_len = if periodic { nbp - 1 } else { nbp };
        let mut rad = vec![0.0f64; rad_len];
        for i in 1..=nbp {
            let coord = par_and_rad[l + i - 1];
            let (x, y) = (coord.x, coord.y);
            par[i - 1] = (uf * (x - wd) + ud * (wf - x)) / (wf - wd);
            if !periodic || i != nbp {
                rad[i - 1] = y;
            }
        }
        let mut inter = LawInterpolate::with_parameters(rad, par, periodic, CONFUSION);
        inter.perform();
        let curve = inter.curve();
        super::law_s::set_curve(&mut self.base, curve);
    }

    /// OCCT Law_Interpol::Set(ParAndRad, Dd, Df, Periodic)
    /// (Law_Interpol.cxx L109-142).
    pub fn set_with_tangents(
        &mut self,
        par_and_rad: &[DVec2],
        dd: f64,
        df: f64,
        periodic: bool,
    ) {
        let l = 0usize; // OCCT: ParAndRad.Lower()
        let nbp = par_and_rad.len(); // OCCT: ParAndRad.Length()

        let mut par = vec![0.0f64; nbp];
        let rad_len = if periodic { nbp - 1 } else { nbp };
        let mut rad = vec![0.0f64; rad_len];
        for i in 1..=nbp {
            let coord = par_and_rad[l + i - 1];
            let (x, y) = (coord.x, coord.y);
            par[i - 1] = x;
            if !periodic || i != nbp {
                rad[i - 1] = y;
            }
        }
        let mut inter = LawInterpolate::with_parameters(rad, par, periodic, CONFUSION);
        // OCCT L139: inter.Load(Dd, Df);
        inter.load_end_tangents(dd, df);
        inter.perform();
        let curve = inter.curve();
        super::law_s::set_curve(&mut self.base, curve);
    }

    /// OCCT Law_Interpol::SetInRelative(ParAndRad, Ud, Uf, Dd, Df, Periodic)
    /// (Law_Interpol.cxx L146-181).
    pub fn set_in_relative_with_tangents(
        &mut self,
        par_and_rad: &[DVec2],
        ud: f64,
        uf: f64,
        dd: f64,
        df: f64,
        periodic: bool,
    ) {
        let l = 0usize; // OCCT: ParAndRad.Lower()
        let u = par_and_rad.len() - 1; // OCCT: ParAndRad.Upper()
        let wd = par_and_rad[l].x;
        let wf = par_and_rad[u].x;
        let nbp = par_and_rad.len(); // OCCT: ParAndRad.Length()

        let mut par = vec![0.0f64; nbp];
        let rad_len = if periodic { nbp - 1 } else { nbp };
        let mut rad = vec![0.0f64; rad_len];
        for i in 1..=nbp {
            let coord = par_and_rad[l + i - 1];
            let (x, y) = (coord.x, coord.y);
            par[i - 1] = (uf * (x - wd) + ud * (wf - x)) / (wf - wd);
            if !periodic || i != nbp {
                rad[i - 1] = y;
            }
        }
        let mut inter = LawInterpolate::with_parameters(rad, par, periodic, CONFUSION);
        // OCCT L178: inter.Load(Dd, Df);
        inter.load_end_tangents(dd, df);
        inter.perform();
        let curve = inter.curve();
        super::law_s::set_curve(&mut self.base, curve);
    }
}

// OCCT inheritance: every Law_Function virtual of Law_BSpFunc applies to
// Law_Interpol; the impls below are pure delegation to the base subobject.
impl LawFunction for LawInterpol {
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
