//! OCCT HelixGeom_BuilderHelix (TKHelix HelixGeom package).
//!
//! 1:1 translation of `HelixGeom_BuilderHelix.hxx` + `.cxx` (L24-172) —
//! builds helix curves using an arbitrary axis, segmenting by full turns and
//! reusing the first coil for cylindrical helixes.

use super::builder_helix_coil::BuilderHelixCoil;
use super::builder_helix_gen::BuilderHelixGen;
use rcad_kernel::geom::BSplineCurve3;
use rcad_kernel::math::gp::{Ax2, Ax3, Trsf};

/// OCCT HelixGeom_BuilderHelix.
#[derive(Debug, Clone)]
pub struct BuilderHelix {
    /// OCCT inheritance chain: BuilderApproxCurve -> BuilderHelixGen.
    pub hgen: BuilderHelixGen,
    pub my_position: Ax2,
}

impl BuilderHelix {
    /// OCCT HelixGeom_BuilderHelix() (BuilderHelix.cxx L24-26).
    pub fn new() -> Self {
        BuilderHelix {
            hgen: BuilderHelixGen::new(),
            my_position: Ax2::new(
                glam::DVec3::ZERO,
                glam::DVec3::Z,
                glam::DVec3::X,
            ),
        }
    }

    /// OCCT SetPosition (L34-37).
    pub fn set_position(&mut self, a_ax2: &Ax2) {
        self.my_position = *a_ax2;
    }

    /// OCCT Position (L40-43).
    pub fn position(&self) -> &Ax2 {
        &self.my_position
    }

    /// OCCT HelixGeom_BuilderHelix::Perform (L48-172).
    pub fn perform(&mut self) {
        self.hgen.base.my_error_status = 0;
        self.hgen.base.my_warning_status = 0;

        // Initialize result containers.
        self.hgen.base.my_curves.clear();
        self.hgen.base.my_tol_reached = -1.0;
        let a_two_pi = 2.0 * std::f64::consts::PI;

        let mut a_bhc = BuilderHelixCoil::new();
        a_bhc.hgen.base.my_tolerance = self.hgen.base.my_tolerance;
        a_bhc
            .hgen
            .base
            .set_approx_parameters(
                self.hgen.base.my_cont,
                self.hgen.base.my_max_degree,
                self.hgen.base.my_max_seg,
            );

        // Determine number of full turns for segmentation.
        let d_t = self.hgen.my_t2 - self.hgen.my_t1;
        let a_n = (d_t / a_two_pi) as i32;
        if a_n == 0 {
            a_bhc.set_curve_parameters(
                self.hgen.my_t1,
                self.hgen.my_t2,
                self.hgen.my_pitch,
                self.hgen.my_r_start,
                self.hgen.my_taper_angle,
                self.hgen.my_is_clock_wise,
            );
            a_bhc.perform();
            let i_err = a_bhc.error_status();
            if i_err != 0 {
                self.hgen.base.my_error_status = 2;
                return;
            }
            let a_c = a_bhc.curves()[0].clone();
            self.hgen.base.my_curves.push(a_c);
            self.hgen.base.my_tol_reached = a_bhc.tolerance_reached();
        } else {
            // Case: helix spans multiple full turns - process in segments.
            let a_tol_angle = 1.0e-4;
            let b_is_cylindrical = self.hgen.my_taper_angle.abs() < a_tol_angle;
            let mut a_t1x = self.hgen.my_t1;
            let mut a_t2x = self.hgen.my_t1 + a_two_pi;
            for i in 1..=a_n {
                if i > 1 && b_is_cylindrical {
                    // Optimization: for cylindrical helixes, reuse first coil
                    // with translation.
                    let a_c1 = self.hgen.base.my_curves[0].clone();
                    let a_p1 = rcad_kernel::math::bspl::de_boor(
                        a_c1.degree,
                        &a_c1.knots,
                        &a_c1.control_points,
                        &a_c1.weights,
                        a_c1.first_parameter(),
                    );
                    let a_pi = glam::DVec3::new(
                        a_p1.x,
                        a_p1.y,
                        a_p1.z + (i as f64 - 1.0) * self.hgen.my_pitch,
                    );
                    let a_ci = a_c1.translated(a_p1, a_pi);
                    self.hgen.base.my_curves.push(a_ci);
                    a_t1x = a_t2x;
                    a_t2x = a_t1x + a_two_pi;
                    // Skip to next iteration for optimization.
                    continue;
                }

                a_bhc.set_curve_parameters(
                    a_t1x,
                    a_t2x,
                    self.hgen.my_pitch,
                    self.hgen.my_r_start,
                    self.hgen.my_taper_angle,
                    self.hgen.my_is_clock_wise,
                );
                // Perform approximation for this segment.
                a_bhc.perform();
                let i_err = a_bhc.error_status();
                if i_err != 0 {
                    self.hgen.base.my_error_status = 2;
                    return;
                }
                // Extract approximated curves from builder.
                let a_c = a_bhc.curves()[0].clone();
                self.hgen.base.my_curves.push(a_c);
                let a_tr = a_bhc.tolerance_reached();
                if a_tr > self.hgen.base.my_tol_reached {
                    self.hgen.base.my_tol_reached = a_tr;
                }
                // Move to next segment parameters.
                a_t1x = a_t2x;
                a_t2x = a_t1x + a_two_pi;
            } // for (i=1; i<=aN; ++i)
            // Handle remaining partial turn if any.
            a_t2x = self.hgen.my_t2;
            let eps = 1.0e-7 * a_two_pi;
            if (a_t2x - a_t1x).abs() > eps {
                a_bhc.set_curve_parameters(
                    a_t1x,
                    a_t2x,
                    self.hgen.my_pitch,
                    self.hgen.my_r_start,
                    self.hgen.my_taper_angle,
                    self.hgen.my_is_clock_wise,
                );
                a_bhc.perform();
                let i_err = a_bhc.error_status();
                if i_err != 0 {
                    self.hgen.base.my_error_status = 2;
                    return;
                }
                // Extract curves from the final partial segment.
                let a_c = a_bhc.curves()[0].clone();
                self.hgen.base.my_curves.push(a_c);
                let a_tr = a_bhc.tolerance_reached();
                if a_tr > self.hgen.base.my_tol_reached {
                    self.hgen.base.my_tol_reached = a_tr;
                }
            }
        }
        // Apply coordinate system transformation to all curves.
        let a_ax3 = Ax3::new();
        let a_ax3x = Ax3::from_ax2(&self.my_position);
        let a_trsf = Trsf::set_displacement(&a_ax3, &a_ax3x);
        // Apply transformation to all generated curves.
        for a_c in self.hgen.base.my_curves.iter_mut() {
            a_c.transform_trsf(&a_trsf);
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

    /// OCCT BuilderApproxCurve::ErrorStatus (forwarded).
    pub fn error_status(&self) -> i32 {
        self.hgen.base.error_status()
    }

    /// OCCT BuilderApproxCurve::WarningStatus (forwarded).
    pub fn warning_status(&self) -> i32 {
        self.hgen.base.warning_status()
    }

    /// OCCT BuilderApproxCurve::ToleranceReached (forwarded).
    pub fn tolerance_reached(&self) -> f64 {
        self.hgen.base.tolerance_reached()
    }

    /// OCCT BuilderApproxCurve::Curves (forwarded).
    pub fn curves(&self) -> &Vec<BSplineCurve3> {
        self.hgen.base.curves()
    }

    /// OCCT BuilderApproxCurve::SetTolerance (forwarded).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.hgen.base.set_tolerance(tol);
    }

    /// OCCT BuilderApproxCurve::SetApproxParameters (forwarded).
    pub fn set_approx_parameters(&mut self, cont: rcad_kernel::math::GeomAbsShape, max_degree: i32, max_seg: i32) {
        self.hgen.base.set_approx_parameters(cont, max_degree, max_seg);
    }
}
