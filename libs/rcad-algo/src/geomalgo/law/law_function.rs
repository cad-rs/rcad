//! OCCT Law_Function (TKGeomAlgo/Law) — 1:1 port of Law_Function.hxx
//! (root class for evolution laws) + Law_Function.cxx (whole file, empty
//! besides RTTI).
//!
//! Architecture mapping: the `Standard_Transient` inheritance hierarchy is
//! expressed as a Rust trait; OCCT `handle<Law_Function>` (a shared,
//! mutable-through-handle reference) maps to
//! `Rc<RefCell<dyn LawFunction>>`, because OCCT methods like `Value` and
//! `Prepare` are declared non-const and mutate through shared handles.

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::GeomAbsShape;

/// OCCT `handle<Law_Function>`.
pub type LawFunctionHandle = Rc<RefCell<dyn LawFunction>>;

/// OCCT Law_Function — root class for evolution laws
/// (Law_Function.hxx L31-63).
pub trait LawFunction {
    /// OCCT Continuity() — pure virtual.
    fn continuity(&self) -> GeomAbsShape;

    /// OCCT NbIntervals(S) — pure virtual: the number of intervals for
    /// continuity S.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize;

    /// OCCT Intervals(T, S) — pure virtual: the parameters bounding the
    /// intervals of continuity S.
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape);

    /// OCCT Value(X) — pure virtual: the value of the function at X.
    fn value(&mut self, x: f64) -> f64;

    /// OCCT D1(X, F, D) — pure virtual: the value and first derivative.
    fn d1(&mut self, x: f64, f: &mut f64, d: &mut f64);

    /// OCCT D2(X, F, D, D2) — pure virtual: value, first and second
    /// derivatives.
    fn d2(&mut self, x: f64, f: &mut f64, d: &mut f64, d2: &mut f64);

    /// OCCT Trim(PFirst, PLast, Tol) — pure virtual: a law equivalent of
    /// this one between PFirst and PLast.
    fn trim(&self, pfirst: f64, plast: f64, tol: f64) -> LawFunctionHandle;

    /// OCCT Bounds(PFirst, PLast) — pure virtual: the parametric bounds.
    fn bounds(&self, pfirst: &mut f64, plast: &mut f64);
}
