//! GAP carriers for Bisector's external OCCT dependencies that are not yet
//! translated in rcad.  Each carrier keeps the exact OCCT call shape and the
//! OCCT failure-path structure; the owning package translation will replace
//! them.  Not part of the Bisector module proper.
//!
//! Dependencies covered:
//! - TKGeomAlgo/GccInt  (GccInt_Bisec/BLine/BCirc/BPar/BHpr/BEll)
//! - TKGeomAlgo/GccAna  (Circ2dBisec / CircLin2dBisec / CircPnt2dBisec /
//!                        LinPnt2dBisec / Pnt2dBisec)
//! - TKMath             (math_BissecNewton / math_FunctionRoot)
//! - TKGeomBase/Geom2dAPI (Geom2dAPI_ProjectPointOnCurve)

use std::cell::RefCell;

use glam::DVec2;
use rcad_kernel::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};
use rcad_kernel::math::root::FunctionWithDerivative;
use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};

use crate::geomalgo::extrema_gen_ext_pc2d::EPCOfExtPC2d;
use crate::geomalgo::geom2d_int::Curve2dAdaptor;

// ===========================================================================
// TKGeomAlgo/GccInt — GccInt_Bisec + the typed results.
// OCCT: GccInt_Bisec.hxx, GccInt_BLine.hxx, GccInt_BCirc.hxx,
//       GccInt_BParab.hxx, GccInt_BHprb.hxx, GccInt_BElips.hxx,
//       GccInt_IType.hxx.
// ===========================================================================

/// OCCT GccInt_IType (GccInt_IType.hxx) — GccInt_Lin/Cir/Ell/Par/Hpr.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GccIntIType {
    Lin,
    Cir,
    Ell,
    Par,
    Hpr,
}

/// OCCT GccInt_Bisec — the typed bisector result.  The OCCT sub-classes
/// (BLine/BCirc/BParab/BHprb/BElips) map to the payload variants.
#[derive(Debug, Clone, Copy)]
pub enum GccIntBisec {
    /// OCCT GccInt_BLine — the gp_Lin2d payload.
    Line(Line2d),
    /// OCCT GccInt_BCirc — the gp_Circ2d payload.
    Circle(Circle2d),
    /// OCCT GccInt_BParab — the gp_Parab2d payload.
    Parabola(Parabola2d),
    /// OCCT GccInt_BHprb — the gp_Hypr2d payload.
    Hyperbola(Hyperbola2d),
    /// OCCT GccInt_BElips — the gp_Elips2d payload.
    Ellipse(Ellipse2d),
}

impl GccIntBisec {
    /// OCCT GccInt_Bisec::ArcType().
    pub fn arc_type(&self) -> GccIntIType {
        match self {
            GccIntBisec::Line(_) => GccIntIType::Lin,
            GccIntBisec::Circle(_) => GccIntIType::Cir,
            GccIntBisec::Parabola(_) => GccIntIType::Par,
            GccIntBisec::Hyperbola(_) => GccIntIType::Hpr,
            GccIntBisec::Ellipse(_) => GccIntIType::Ell,
        }
    }

    /// OCCT GccInt_BLine::Line() (raises on the other types).
    pub fn line(&self) -> Line2d {
        match self {
            GccIntBisec::Line(l) => *l,
            _ => unreachable!("GccInt_Bisec::Line() type mismatch"),
        }
    }

    /// OCCT GccInt_BCirc::Circle().
    pub fn circle(&self) -> Circle2d {
        match self {
            GccIntBisec::Circle(c) => *c,
            _ => unreachable!("GccInt_Bisec::Circle() type mismatch"),
        }
    }

    /// OCCT GccInt_BParab::Parabola().
    pub fn parabola(&self) -> Parabola2d {
        match self {
            GccIntBisec::Parabola(p) => *p,
            _ => unreachable!("GccInt_Bisec::Parabola() type mismatch"),
        }
    }

    /// OCCT GccInt_BHprb::Hyperbola().
    pub fn hyperbola(&self) -> Hyperbola2d {
        match self {
            GccIntBisec::Hyperbola(h) => *h,
            _ => unreachable!("GccInt_Bisec::Hyperbola() type mismatch"),
        }
    }

    /// OCCT GccInt_BElips::Ellipse().
    pub fn ellipse(&self) -> Ellipse2d {
        match self {
            GccIntBisec::Ellipse(e) => *e,
            _ => unreachable!("GccInt_Bisec::Ellipse() type mismatch"),
        }
    }
}

// ===========================================================================
// TKGeomAlgo/GccAna — the five analytic bisector constructors used by
// Bisector_BisecAna.  Carriers preserve the OCCT interface and the OCCT
// failure path: not-done / no solutions -> BisecAna::Perform leaves
// `thebisector` unset (OCCT: null handle).
// ===========================================================================

/// GAP: OCCT GccAna_Circ2dBisec (GccAna_Circ2dBisec.cxx) — circle/circle
/// bisector (lines + circles).  Not yet translated; carries the not-done
/// failure path.
pub struct GccAnaCirc2dBisec {
    _circle1: Circle2d,
    _circle2: Circle2d,
}

impl GccAnaCirc2dBisec {
    /// OCCT GccAna_Circ2dBisec(Circle1, Circle2, Tolerance).
    pub fn new(circle1: Circle2d, circle2: Circle2d, _tolerance: f64) -> Self {
        GccAnaCirc2dBisec { _circle1: circle1, _circle2: circle2 }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT NbSolutions().
    pub fn nb_solutions(&self) -> usize {
        0
    }

    /// OCCT ThisSolution(Index).
    pub fn this_solution(&self, _index: usize) -> GccIntBisec {
        panic!("GAP: GccAna_Circ2dBisec not yet translated")
    }
}

/// GAP: OCCT GccAna_CircLin2dBisec (GccAna_CircLin2dBisec.cxx) —
/// circle/line bisector (line + parabolas).  Not yet translated.
pub struct GccAnaCircLin2dBisec {
    _circle: Circle2d,
    _line: Line2d,
}

impl GccAnaCircLin2dBisec {
    /// OCCT GccAna_CircLin2dBisec(Circle, Line).
    pub fn new(circle: Circle2d, line: Line2d) -> Self {
        GccAnaCircLin2dBisec { _circle: circle, _line: line }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT NbSolutions().
    pub fn nb_solutions(&self) -> usize {
        0
    }

    /// OCCT ThisSolution(Index).
    pub fn this_solution(&self, _index: usize) -> GccIntBisec {
        panic!("GAP: GccAna_CircLin2dBisec not yet translated")
    }
}

/// GAP: OCCT GccAna_CircPnt2dBisec (GccAna_CircPnt2dBisec.cxx) —
/// circle/point bisector (circles, parabolas, hyperbolas, ellipse, line).
/// Not yet translated.
pub struct GccAnaCircPnt2dBisec {
    _circle: Circle2d,
    _point: DVec2,
    _tolerance: f64,
}

impl GccAnaCircPnt2dBisec {
    /// OCCT GccAna_CircPnt2dBisec(Circle, Point, Tolerance).
    pub fn new(circle: Circle2d, point: DVec2, tolerance: f64) -> Self {
        GccAnaCircPnt2dBisec { _circle: circle, _point: point, _tolerance: tolerance }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT NbSolutions().
    pub fn nb_solutions(&self) -> usize {
        0
    }

    /// OCCT ThisSolution(Index).
    pub fn this_solution(&self, _index: usize) -> GccIntBisec {
        panic!("GAP: GccAna_CircPnt2dBisec not yet translated")
    }
}

/// GAP: OCCT GccAna_LinPnt2dBisec (GccAna_LinPnt2dBisec.cxx) — line/point
/// bisector (a parabola, or a line when the point is on the line).
/// Not yet translated; OCCT IsDone() is true after construction — the
/// carrier preserves the OCCT raise path of ThisSolution().
pub struct GccAnaLinPnt2dBisec {
    _line: Line2d,
    _point: DVec2,
}

impl GccAnaLinPnt2dBisec {
    /// OCCT GccAna_LinPnt2dBisec(Line, Point).
    pub fn new(line: Line2d, point: DVec2) -> Self {
        GccAnaLinPnt2dBisec { _line: line, _point: point }
    }

    /// OCCT ThisSolution() — the unique bisector.
    pub fn this_solution(&self) -> GccIntBisec {
        panic!("GAP: GccAna_LinPnt2dBisec not yet translated")
    }
}

/// GAP: OCCT GccAna_Pnt2dBisec (GccAna_Pnt2dBisec.cxx) — point/point
/// bisector (the perpendicular line).  Not yet translated; same raise-path
/// preservation as LinPnt2dBisec.
pub struct GccAnaPnt2dBisec {
    _point1: DVec2,
    _point2: DVec2,
}

impl GccAnaPnt2dBisec {
    /// OCCT GccAna_Pnt2dBisec(Point1, Point2).
    pub fn new(point1: DVec2, point2: DVec2) -> Self {
        GccAnaPnt2dBisec { _point1: point1, _point2: point2 }
    }

    /// OCCT ThisSolution() — the perpendicular bisector line.
    pub fn this_solution(&self) -> Line2d {
        panic!("GAP: GccAna_Pnt2dBisec not yet translated")
    }
}

// ===========================================================================
// TKMath — math_BissecNewton / math_FunctionRoot carriers.
// ===========================================================================

/// OCCT math_BissecNewton — bisection+Newton root finder on
/// math_FunctionWithDerivative.  GAP: TKMath math_BissecNewton is not yet
/// translated in rcad; this carrier delegates to the kernel
/// `math::root::biss_newton` port of the same algorithm.
pub struct MathBissecNewton {
    /// OCCT Tol.
    tol: f64,
    /// OCCT Done.
    done: bool,
    /// OCCT Sol.
    sol: f64,
}

impl MathBissecNewton {
    /// OCCT math_BissecNewton(Tolerance).
    pub fn new(tol: f64) -> Self {
        MathBissecNewton { tol, done: false, sol: 0.0 }
    }

    /// OCCT Perform(F, Bounds1, Bounds2, NbIterations).
    pub fn perform(
        &mut self,
        f: &mut dyn FunctionWithDerivative,
        bound1: f64,
        bound2: f64,
        nb_iterations: i32,
    ) {
        // Shim: the kernel biss_newton implements the bissec+Newton search;
        // NbIterations maps to its fixed iteration budget.
        let _ = nb_iterations;
        let fcell = RefCell::new(f);
        let value = |x: f64| fcell.borrow_mut().value(x).unwrap_or(f64::MAX);
        let deriv = |x: f64| fcell.borrow_mut().derivative(x).unwrap_or(1.0);
        if let Some(sol) = rcad_kernel::math::root::biss_newton(value, deriv, bound1, bound2, self.tol) {
            self.done = true;
            self.sol = sol;
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Root().
    pub fn root(&self) -> f64 {
        self.sol
    }
}

/// OCCT math_FunctionRoot — Newton iteration from a guess, clamped to
/// [Binf, Bsup].  GAP: TKMath math_FunctionRoot is not yet translated in
/// rcad; the carrier reproduces the OCCT Perform loop
/// (math_FunctionRoot.cxx — Newton steps with bound clamping, convergence
/// on EpsX/EpsF).
pub struct MathFunctionRoot {
    done: bool,
    root: f64,
}

impl MathFunctionRoot {
    /// OCCT math_FunctionRoot(F, Guess, EpsX, Binf, Bsup) — the bounded
    /// constructor; NbIterations defaults to 100 in OCCT.
    pub fn new(
        f: &mut dyn FunctionWithDerivative,
        guess: f64,
        eps_x: f64,
        binf: f64,
        bsup: f64,
    ) -> Self {
        Self::new_with_iterations(f, guess, eps_x, binf, bsup, 100)
    }

    /// OCCT math_FunctionRoot(F, Guess, EpsX, Binf, Bsup, NbIterations).
    pub fn new_with_iterations(
        f: &mut dyn FunctionWithDerivative,
        guess: f64,
        eps_x: f64,
        binf: f64,
        bsup: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut done = false;
        let mut root = guess;
        let (lo, hi) = if binf <= bsup { (binf, bsup) } else { (bsup, binf) };
        let mut x = guess.clamp(lo, hi);
        for _ in 0..nb_iterations.max(1) {
            let (fv, dv) = match f.values(x) {
                Some((v, d)) => (v, d),
                None => break,
            };
            root = x;
            if fv.abs() < eps_x || dv.abs() < f64::MIN_POSITIVE {
                if fv.abs() < eps_x {
                    done = true;
                }
                break;
            }
            let next = x - fv / dv;
            if (next - x).abs() < eps_x {
                x = next.clamp(lo, hi);
                root = x;
                done = true;
                break;
            }
            x = next.clamp(lo, hi);
        }
        MathFunctionRoot { done, root }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Root().
    pub fn root(&self) -> f64 {
        self.root
    }
}

// ===========================================================================
// TKGeomBase/Geom2dAPI — Geom2dAPI_ProjectPointOnCurve (the 2D API).
// OCCT engine: Extrema_ExtPC2d(P, C, U1, U2) -> lower-distance point.
// GAP: Geom2dAPI is not translated; this carrier drives the rcad
// Extrema_ExtPC2d port (EPCOfExtPC2d) directly.
// ===========================================================================

pub struct Geom2dAPIProjectPointOnCurve {
    /// The lower-distance extremum parameter (OCCT: theNbSolutions points,
    /// of which only the closest is used by Bisector).
    nb_points: usize,
    lower_distance_parameter: f64,
}

impl Geom2dAPIProjectPointOnCurve {
    /// OCCT Geom2dAPI_ProjectPointOnCurve(P, C, U1, U2).
    pub fn new(p: DVec2, c: &dyn Curve2dAdaptor, u1: f64, u2: f64) -> Self {
        // OCCT Extrema_ExtPC2d(P, C, NbSamples=20?, TolU, TolF) defaults.
        let mut ext = EPCOfExtPC2d::with_range(p, c, 20, u1, u2, PCONFUSION, CONFUSION);
        ext.perform(p);
        let mut nb_points = 0usize;
        let mut best_param = u1;
        let mut best_dist2 = f64::INFINITY;
        if ext.is_done() {
            for n in 0..ext.nb_ext() {
                if !ext.is_min(n) {
                    continue;
                }
                nb_points += 1;
                let d2 = ext.square_distance(n);
                if d2 < best_dist2 {
                    best_dist2 = d2;
                    best_param = ext.point(n).parameter();
                }
            }
        }
        Geom2dAPIProjectPointOnCurve {
            nb_points,
            lower_distance_parameter: best_param,
        }
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.nb_points
    }

    /// OCCT LowerDistanceParameter().
    pub fn lower_distance_parameter(&self) -> f64 {
        self.lower_distance_parameter
    }
}
