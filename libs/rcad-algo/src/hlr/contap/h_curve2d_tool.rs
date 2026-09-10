// OCCT Contap_HCurve2dTool (TKHLR) — the static tool over the contour 2D
// arc (the `TheArcTool` of Contap_TheSearch).
//
// Contap_HCurve2dTool.hxx L35-131 + Contap_HCurve2dTool.lxx L34-195 (each
// one-liner delegates to the Adaptor2d_Curve2d) + Contap_HCurve2dTool.cxx
// L21-52 (NbSamples).  The OCCT `occ::handle<Adaptor2d_Curve2d>` maps to
// `&dyn Curve2dAdaptor`.

use glam::DVec2;

use crate::geomalgo::geom2d_int::{Curve2dAdaptor, Curve2dType};

/// The 2D arc type used across Contap — the OCCT
/// `occ::handle<Adaptor2d_Curve2d>`.
pub type HCurve2d = std::sync::Arc<dyn Curve2dAdaptor>;

/// OCCT Contap_HCurve2dTool::FirstParameter (lxx L35-38).
pub fn first_parameter(c: &dyn Curve2dAdaptor) -> f64 {
    c.first_parameter()
}

/// OCCT Contap_HCurve2dTool::LastParameter (lxx L41-44).
pub fn last_parameter(c: &dyn Curve2dAdaptor) -> f64 {
    c.last_parameter()
}

/// OCCT Contap_HCurve2dTool::IsClosed (lxx L68-71).
pub fn is_closed(c: &dyn Curve2dAdaptor) -> bool {
    c.is_closed()
}

/// OCCT Contap_HCurve2dTool::IsPeriodic (lxx L74-77).
pub fn is_periodic(c: &dyn Curve2dAdaptor) -> bool {
    c.is_periodic()
}

/// OCCT Contap_HCurve2dTool::Period (lxx L80-83).
pub fn period(c: &dyn Curve2dAdaptor) -> f64 {
    c.period()
}

/// OCCT Contap_HCurve2dTool::Value (lxx L86-89).
pub fn value(c: &dyn Curve2dAdaptor, u: f64) -> DVec2 {
    c.value(u)
}

/// OCCT Contap_HCurve2dTool::D0 (lxx L92-97).
pub fn d0(c: &dyn Curve2dAdaptor, u: f64) -> DVec2 {
    c.value(u)
}

/// OCCT Contap_HCurve2dTool::D1 (lxx L100-106).
pub fn d1(c: &dyn Curve2dAdaptor, u: f64) -> (DVec2, DVec2) {
    c.d1(u)
}

/// OCCT Contap_HCurve2dTool::D2 (lxx L109-117).
pub fn d2(c: &dyn Curve2dAdaptor, u: f64) -> (DVec2, DVec2, DVec2) {
    c.d2(u)
}

/// OCCT Contap_HCurve2dTool::DN (lxx L132-138).
pub fn dn(c: &dyn Curve2dAdaptor, u: f64, n: i32) -> DVec2 {
    c.dn(u, n)
}

/// OCCT Contap_HCurve2dTool::Resolution (lxx L141-145).
pub fn resolution(c: &dyn Curve2dAdaptor, r3d: f64) -> f64 {
    c.resolution(r3d)
}

/// OCCT Contap_HCurve2dTool::GetType (lxx L148-151).
pub fn get_type(c: &dyn Curve2dAdaptor) -> Curve2dType {
    c.get_type()
}

/// OCCT Contap_HCurve2dTool::NbSamples (cxx L21-52).
pub fn nb_samples(c: &dyn Curve2dAdaptor, u0: f64, u1: f64) -> i32 {
    let mut nbs = 10.0;
    match c.get_type() {
        Curve2dType::Line => {
            nbs = 2.0;
        }
        Curve2dType::BezierCurve => {
            nbs = 3.0 + c.nb_poles() as f64;
        }
        Curve2dType::BSplineCurve => {
            nbs = c.nb_knots() as f64;
            nbs *= c.degree() as f64;
            nbs *= c.last_parameter() - c.first_parameter();
            nbs /= u1 - u0;
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        _ => {}
    }
    if nbs > 50.0 {
        nbs = 50.0;
    }
    nbs as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle2d, Curve2d};

    /// OCCT anchor: NbSamples dispatch — Line gives 2, Circle keeps the
    /// default 10, and the result is clamped at 50 (cxx L21-52).
    #[test]
    fn hcurve2d_tool_nb_samples() {
        let line = Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::new(1.0, 0.0),
        });
        assert_eq!(nb_samples(&line, 0.0, 1.0), 2);

        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: 1.0,
        });
        assert_eq!(nb_samples(&circle, 0.0, 1.0), 10);
        // Value/D1 round-trip on the circle.
        let (p, t) = d1(&circle, 0.0);
        assert!((p.y - 0.0).abs() < 1e-12);
        assert!(t.x.abs() < 1e-9 && (t.y - 1.0).abs() < 1e-9);
        // NB: the enum-level Curve2dEval::is_periodic is not delegated to
        // the variants in the kernel (returns the trait default false); the
        // Curve2dAdaptor period() carries the OCCT semantics.
        assert!((period(&circle) - std::f64::consts::TAU).abs() < 1e-12);
    }
}
