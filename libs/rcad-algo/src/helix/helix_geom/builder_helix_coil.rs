//! OCCT HelixGeom_BuilderHelixCoil (TKHelix HelixGeom package).
//!
//! 1:1 translation of `HelixGeom_BuilderHelixCoil.hxx` + `.cxx` — builds one
//! helix coil (a single 2π-range BSpline approximation) with axis OZ.

use super::builder_helix_gen::BuilderHelixGen;
use super::tools;
use rcad_kernel::geom::BSplineCurve3;

/// OCCT HelixGeom_BuilderHelixCoil.
#[derive(Debug, Clone, Default)]
pub struct BuilderHelixCoil {
    /// OCCT inheritance chain: BuilderApproxCurve -> BuilderHelixGen.
    pub hgen: BuilderHelixGen,
}

impl BuilderHelixCoil {
    /// OCCT HelixGeom_BuilderHelixCoil() (BuilderHelixCoil.cxx L21-30).
    pub fn new() -> Self {
        BuilderHelixCoil {
            hgen: BuilderHelixGen::new(),
        }
    }

    /// OCCT HelixGeom_BuilderHelixCoil::Perform
    /// (HelixGeom_BuilderHelixCoil.cxx L38-68).
    pub fn perform(&mut self) {
        self.hgen.base.my_error_status = 0;
        self.hgen.base.my_warning_status = 0;
        // Initialize variables for curve approximation.
        // Clear previous results and setup helix adaptor.
        self.hgen.base.my_curves.clear();
        // Load helix parameters into the adaptor.
        let mut a_adaptor = super::helix_curve::HelixCurve::new();
        a_adaptor.load(
            self.hgen.my_t1,
            self.hgen.my_t2,
            self.hgen.my_pitch,
            self.hgen.my_r_start,
            self.hgen.my_taper_angle,
            self.hgen.my_is_clock_wise,
        );
        // Perform B-spline approximation of the helix curve (OCCT L91-97):
        // myTolReached is passed as the ApprCurve3D output reference.
        let (i_err, a_bc, the_max_error) = tools::appr_curve3d(
            &a_adaptor,
            self.hgen.base.my_tolerance,
            self.hgen.base.my_cont,
            self.hgen.base.my_max_seg,
            self.hgen.base.my_max_degree,
        );
        // OCCT ApprCurve3D writes theMaxError from its L139 onward, i.e. on
        // the return-2 and return-0 paths; the return-1 path (approximation
        // not done) leaves it untouched.
        if i_err != 1 {
            self.hgen.base.my_tol_reached = the_max_error;
        }
        if i_err != 0 {
            self.hgen.base.my_error_status = 2;
        } else {
            self.hgen
                .base
                .my_curves
                .push(a_bc.expect("ApprCurve3D returned 0 without a BSpline"));
        }
    }

    /// OCCT BuilderHelixGen::SetCurveParameters (forwarded).
    pub fn set_curve_parameters(
        &mut self,
        a_t1: f64,
        a_t2: f64,
        a_pitch: f64,
        a_r_start: f64,
        a_taper_angle: f64,
        a_is_cw: bool,
    ) {
        self.hgen
            .set_curve_parameters(a_t1, a_t2, a_pitch, a_r_start, a_taper_angle, a_is_cw);
    }

    /// OCCT BuilderHelixGen::CurveParameters (forwarded).
    #[allow(clippy::type_complexity)]
    pub fn curve_parameters(&self) -> (f64, f64, f64, f64, f64, bool) {
        self.hgen.curve_parameters()
    }

    /// OCCT BuilderApproxCurve::SetApproxParameters (forwarded).
    pub fn set_approx_parameters(&mut self, a_cont: rcad_kernel::math::GeomAbsShape, a_max_degree: i32, a_max_seg: i32) {
        self.hgen.base.set_approx_parameters(a_cont, a_max_degree, a_max_seg);
    }

    /// OCCT BuilderApproxCurve::ApproxParameters (forwarded).
    pub fn approx_parameters(&self) -> (rcad_kernel::math::GeomAbsShape, i32, i32) {
        self.hgen.base.approx_parameters()
    }

    /// OCCT BuilderApproxCurve::SetTolerance (forwarded).
    pub fn set_tolerance(&mut self, a_tolerance: f64) {
        self.hgen.base.set_tolerance(a_tolerance);
    }

    /// OCCT BuilderApproxCurve::Tolerance (forwarded).
    pub fn tolerance(&self) -> f64 {
        self.hgen.base.tolerance()
    }

    /// OCCT BuilderApproxCurve::ErrorStatus (forwarded).
    pub fn error_status(&self) -> i32 {
        self.hgen.base.error_status()
    }

    /// OCCT BuilderApproxCurve::ToleranceReached (forwarded).
    pub fn tolerance_reached(&self) -> f64 {
        self.hgen.base.tolerance_reached()
    }

    /// OCCT BuilderApproxCurve::Curves (forwarded).
    pub fn curves(&self) -> &Vec<BSplineCurve3> {
        self.hgen.base.curves()
    }
}
