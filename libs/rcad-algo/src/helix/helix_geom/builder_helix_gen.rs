//! OCCT HelixGeom_BuilderHelixGen (TKHelix HelixGeom package).
//!
//! 1:1 translation of `HelixGeom_BuilderHelixGen.hxx` + `.cxx` — base class
//! adding helix-specific curve parameters over BuilderApproxCurve.

use super::builder_approx_curve::BuilderApproxCurveBase;

/// OCCT HelixGeom_BuilderHelixGen.
#[derive(Debug, Clone, Default)]
pub struct BuilderHelixGen {
    /// OCCT inheritance: HelixGeom_BuilderApproxCurve base part.
    pub base: BuilderApproxCurveBase,
    pub my_t1: f64,
    pub my_t2: f64,
    pub my_pitch: f64,
    pub my_r_start: f64,
    pub my_taper_angle: f64,
    pub my_is_clock_wise: bool,
}

impl BuilderHelixGen {
    /// OCCT HelixGeom_BuilderHelixGen() (BuilderHelixGen.cxx L20-28).
    pub fn new() -> Self {
        BuilderHelixGen {
            base: BuilderApproxCurveBase::new(),
            my_t1: 0.0,
            my_t2: 2.0 * std::f64::consts::PI,
            my_pitch: 1.0,
            my_r_start: 1.0,
            my_taper_angle: 0.0,
            my_is_clock_wise: true,
        }
    }

    /// OCCT SetCurveParameters (L36-51).
    pub fn set_curve_parameters(
        &mut self,
        a_t1: f64,
        a_t2: f64,
        a_pitch: f64,
        a_r_start: f64,
        a_taper_angle: f64,
        a_is_cw: bool,
    ) {
        self.my_t1 = a_t1;
        self.my_t2 = a_t2;
        self.my_pitch = a_pitch;
        self.my_r_start = a_r_start;
        self.my_taper_angle = a_taper_angle;
        self.my_is_clock_wise = a_is_cw;
    }

    /// OCCT CurveParameters (L55-68).
    #[allow(clippy::type_complexity)]
    pub fn curve_parameters(&self) -> (f64, f64, f64, f64, f64, bool) {
        (
            self.my_t1,
            self.my_t2,
            self.my_pitch,
            self.my_r_start,
            self.my_taper_angle,
            self.my_is_clock_wise,
        )
    }
}
