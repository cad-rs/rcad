//! OCCT MAT2d_CutCurve — cuts a curve at the extremas of curvature and at
//! the inflections.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_CutCurve.hxx L17-71, MAT2d_CutCurve.cxx L28-116
//!
//! Note: the .hxx also declares `Perform(C, aSide)` (L49) and `PerformInf(C)`
//! (L52) but the .cxx defines neither (exported declarations without a
//! definition; nothing in OCCT links them).  No body exists to translate.

use rcad_kernel::geom::{Curve2d, TrimmedCurve2};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};

/// GAP carrier (dependency of another package).
///
/// OCCT dependency: `GeomLProp_CurAndInf2d` (TKGeomBase/GeomLProp,
/// GeomLProp_CurAndInf2d.hxx/.cxx) — computes the extrema of curvature and
/// the inflections of a 2d curve.  The GeomLProp package is not translated
/// yet, so this carrier reproduces the public interface consumed by
/// MAT2d_CutCurve::Perform (cxx L43-52: `Sommets.Perform(C)`,
/// `Sommets.IsDone()`, `Sommets.IsEmpty()`, `Sommets.NbPoints()`,
/// `Sommets.Parameter(i)`) with the OCCT *failure path*: `IsDone()` returns
/// false, so the cutting pass leaves `theCurves` empty and the curve reports
/// UnModified — exactly the behavior of OCCT when the CurAndInf computation
/// does not succeed.  When the GeomLProp port lands, its implementation
/// replaces the body of [`CurAndInf2d::perform`].
#[derive(Debug, Default)]
pub struct CurAndInf2d {
    done: bool,
    parameters: Vec<f64>,
}

impl CurAndInf2d {
    /// OCCT GeomLProp_CurAndInf2d() — empty constructor.
    pub fn new() -> Self {
        CurAndInf2d {
            done: false,
            parameters: Vec::new(),
        }
    }

    /// OCCT GeomLProp_CurAndInf2d::Perform(C) — GAP: see type docs.
    pub fn perform(&mut self, _c: &Curve2d) {
        self.done = false;
        self.parameters.clear();
    }

    /// OCCT GeomLProp_CurAndInf2d::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT GeomLProp_CurAndInf2d::IsEmpty().
    pub fn is_empty(&self) -> bool {
        self.parameters.is_empty()
    }

    /// OCCT GeomLProp_CurAndInf2d::NbPoints().
    pub fn nb_points(&self) -> usize {
        self.parameters.len()
    }

    /// OCCT GeomLProp_CurAndInf2d::Parameter(Index) (1-based).
    pub fn parameter(&self, index: usize) -> f64 {
        self.parameters[index - 1]
    }
}

/// OCCT MAT2d_CutCurve (MAT2d_CutCurve.hxx L34-69) — cuts a curve at the
/// extremas of curvature and at the inflections, constructing a trimmed
/// curve for each interval.
#[derive(Debug, Default, Clone)]
pub struct Mat2dCutCurve {
    /// OCCT NCollection_Sequence<occ::handle<Geom2d_Curve>> theCurves
    /// (the appended curves are Geom2d_TrimmedCurve).
    the_curves: Vec<Curve2d>,
}

impl Mat2dCutCurve {
    /// OCCT MAT2d_CutCurve::MAT2d_CutCurve() (cxx L28).
    pub fn new() -> Self {
        Mat2dCutCurve {
            the_curves: Vec::new(),
        }
    }

    /// OCCT MAT2d_CutCurve::MAT2d_CutCurve(C) (cxx L32-35).
    pub fn new_with_curve(c: &Curve2d) -> Self {
        let mut cut = Mat2dCutCurve::new();
        cut.perform(c);
        cut
    }

    /// OCCT MAT2d_CutCurve::Perform(C) (cxx L39-83) — cuts a curve at the
    /// extremas of curvature and at the inflections.
    pub fn perform(&mut self, c: &Curve2d) {
        // OCCT theCurves.Clear().
        self.the_curves.clear();

        // OCCT L43-49: GeomLProp_CurAndInf2d Sommets; ... constexpr double
        // PTol = Precision::PConfusion() * 10; Tol = Confusion() * 10;
        let mut sommets = CurAndInf2d::new();
        let ptol: f64 = PCONFUSION * 10.0;
        let tol: f64 = CONFUSION * 10.0;
        let mut ya_cut = false;
        sommets.perform(c);

        // OCCT L52: if (Sommets.IsDone() && !Sommets.IsEmpty())
        if sommets.is_done() && !sommets.is_empty() {
            let mut uf = c.first_parameter();
            let ul = c.last_parameter();
            let mut pf = c.value(uf);
            let pl = c.value(ul);

            // OCCT L59-76: for (int i = 1; i <= Sommets.NbPoints(); i++)
            for i in 1..=sommets.nb_points() {
                let uc = sommets.parameter(i);

                let pc = c.value(uc);
                if uc - uf > ptol && pc.distance(pf) > tol {
                    if ul - uc < ptol || pl.distance(pc) < tol {
                        break;
                    }
                    // OCCT L70: TrimC = new Geom2d_TrimmedCurve(C, UF, UC).
                    let trim_c = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(c.clone()),
                        t_min: uf,
                        t_max: uc,
                    });
                    self.the_curves.push(trim_c);
                    uf = uc;
                    pf = pc;
                    ya_cut = true;
                }
            }
            // OCCT L77-81.
            if ya_cut {
                let trim_c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c.clone()),
                    t_min: uf,
                    t_max: ul,
                });
                self.the_curves.push(trim_c);
            }
        }
    }

    /// OCCT MAT2d_CutCurve::UnModified() const (cxx L87-90) — true if the
    /// curve is not cut.
    pub fn un_modified(&self) -> bool {
        self.the_curves.is_empty()
    }

    /// OCCT MAT2d_CutCurve::NbCurves() const (cxx L94-101) — raises
    /// Standard_OutOfRange if the curve is UnModified.
    pub fn nb_curves(&self) -> usize {
        if self.un_modified() {
            // OCCT throw Standard_OutOfRange().
            panic!("Standard_OutOfRange: MAT2d_CutCurve::NbCurves on an unmodified curve");
        }
        self.the_curves.len()
    }

    /// OCCT MAT2d_CutCurve::Value(Index) const (cxx L105-116) — returns the
    /// Indexth curve (1-based; raises Standard_OutOfRange out of range).
    pub fn value(&self, index: usize) -> &Curve2d {
        if self.un_modified() {
            panic!("Standard_OutOfRange: MAT2d_CutCurve::Value on an unmodified curve");
        }
        if index < 1 || index > self.the_curves.len() {
            panic!("Standard_OutOfRange: MAT2d_CutCurve::Value index out of range");
        }
        // OCCT: occ::down_cast<Geom2d_TrimmedCurve>(theCurves.Value(Index))
        // — every appended curve is a Geom2d_TrimmedCurve.
        &self.the_curves[index - 1]
    }
}
