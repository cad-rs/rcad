//! OCCT HelixGeom_BuilderApproxCurve (TKHelix HelixGeom package).
//!
//! 1:1 translation of `HelixGeom_BuilderApproxCurve.hxx` + `.cxx` — base
//! class for helix curve approximation algorithms.  C++ inheritance is
//! modeled by embedding `BuilderApproxCurveBase` in the derived builders
//! (`HelixGeom_BuilderHelixGen` -> `HelixGeom_BuilderHelix` /
//! `HelixGeom_BuilderHelixCoil`), with forwarding accessors.

use rcad_kernel::geom::BSplineCurve3;
use rcad_kernel::math::GeomAbsShape;

/// OCCT HelixGeom_BuilderApproxCurve (protected base).
#[derive(Debug, Clone)]
pub(crate) struct BuilderApproxCurveBase {
    pub my_error_status: i32,
    pub my_warning_status: i32,
    pub my_tolerance: f64,
    pub my_cont: GeomAbsShape,
    pub my_max_degree: i32,
    pub my_max_seg: i32,
    pub my_tol_reached: f64,
    /// NCollection_Sequence<handle<Geom_Curve>> — all produced curves are
    /// BSplines.
    pub my_curves: Vec<BSplineCurve3>,
}

impl Default for BuilderApproxCurveBase {
    fn default() -> Self {
        Self::new()
    }
}

impl BuilderApproxCurveBase {
    /// OCCT HelixGeom_BuilderApproxCurve() (BuilderApproxCurve.cxx L20-29).
    pub fn new() -> Self {
        BuilderApproxCurveBase {
            my_error_status: 0,
            my_warning_status: 0,
            my_tolerance: 0.0001,
            my_cont: GeomAbsShape::C2,
            my_max_degree: 8,
            my_max_seg: 150,
            my_tol_reached: 99.0,
            my_curves: Vec::new(),
        }
    }

    /// OCCT SetApproxParameters.
    pub fn set_approx_parameters(&mut self, a_cont: GeomAbsShape, a_max_degree: i32, a_max_seg: i32) {
        self.my_cont = a_cont;
        self.my_max_degree = a_max_degree;
        self.my_max_seg = a_max_seg;
    }

    /// OCCT ApproxParameters.
    pub fn approx_parameters(&self) -> (GeomAbsShape, i32, i32) {
        (self.my_cont, self.my_max_degree, self.my_max_seg)
    }

    /// OCCT SetTolerance.
    pub fn set_tolerance(&mut self, a_tolerance: f64) {
        self.my_tolerance = a_tolerance;
    }

    /// OCCT Tolerance.
    pub fn tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT ToleranceReached.
    pub fn tolerance_reached(&self) -> f64 {
        self.my_tol_reached
    }

    /// OCCT Curves.
    pub fn curves(&self) -> &Vec<BSplineCurve3> {
        &self.my_curves
    }

    /// OCCT ErrorStatus.
    pub fn error_status(&self) -> i32 {
        self.my_error_status
    }

    /// OCCT WarningStatus.
    pub fn warning_status(&self) -> i32 {
        self.my_warning_status
    }
}
