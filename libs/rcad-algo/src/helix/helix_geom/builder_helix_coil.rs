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

    /// OCCT HelixGeom_BuilderHelixCoil::Perform (L38-68).
    pub fn perform(&mut self) {
        self.hgen.base.my_error_status = 0;
        self.hgen.base.my_warning_status = 0;
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
        // Perform B-spline approximation of the helix curve.
        let (i_err, a_bc, tol_reached) = tools::appr_curve3d(
            &a_adaptor,
            self.hgen.base.my_tolerance,
            self.hgen.base.my_cont,
            self.hgen.base.my_max_seg,
            self.hgen.base.my_max_degree,
        );
        if i_err != 0 {
            self.hgen.base.my_error_status = 2;
        } else {
            if let Some(bc) = a_bc {
                self.hgen.base.my_curves.push(bc);
            }
            self.hgen.base.my_tol_reached = tol_reached;
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
