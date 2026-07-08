//! OCCT-aligned TKG3d GTest translations.
//!
//! OCCT source: src/ModelingData/TKG3d/GTests/
//!
//! Tests for 3D curve types: OffsetCurve3, SineWave3, CircularHelix3
//! and surface types: BezierSurface, EllipsoidalSurface, HelicoidSurface
//! and algorithms: extrema_curve_curve, interpolate_points.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;

const TOL: f64 = 1e-7;

// =============================================================================
// Geom_BezierSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod bezier_surface_tests {
    use super::*;

    fn make_bezier_plane() -> Surface3 {
        // A 2x2 control point grid defines a bilinear surface.
        // u in [0,1], v in [0,1]. Poles at corners of 1x1 square in XY.
        Surface3::Bezier(BezierSurface {
            control_points: vec![
                vec![DVec3::ZERO,         DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0,1.0,0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        })
    }

    #[test]
    fn bezier_eval_corner_u0v0() {
        let s = make_bezier_plane();
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn bezier_eval_corner_u1v1() {
        let s = make_bezier_plane();
        let p = s.point_at(1.0, 1.0);
        assert!((p - DVec3::new(1.0, 1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn bezier_eval_midpoint() {
        let s = make_bezier_plane();
        let p = s.point_at(0.5, 0.5);
        assert!((p - DVec3::new(0.5, 0.5, 0.0)).length() < TOL);
    }

    #[test]
    fn bezier_default_domain() {
        let s = make_bezier_plane();
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - 1.0).abs() < TOL);
        assert!((v0 - 0.0).abs() < TOL);
        assert!((v1 - 1.0).abs() < TOL);
    }

    #[test]
    fn bezier_normal_at_mid() {
        let s = make_bezier_plane();
        let n = s.normal_at(0.5, 0.5);
        // Planar surface in XY → normal should be (0,0,±1)
        assert!(n.cross(DVec3::Z).length() < TOL);
    }
}

// =============================================================================
// Geom_OffsetCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod offset_curve_tests {
    use super::*;

    fn make_offset_line() -> Curve3 {
        let base = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        Curve3::Offset(OffsetCurve3 {
            basis: Box::new(base),
            offset_distance: 1.0,
            offset_dir: DVec3::Z,
        })
    }

    #[test]
    fn offset_point_at_zero() {
        let c = make_offset_line();
        // Line along X, tangent = X, offset_dir = Z.
        // offset = offset_distance * (tangent × offset_dir).normalize()
        //        = 1.0 * (X × Z).normalize() = 1.0 * (-Y)
        // point at t=0 = (0,0,0) + (0,-1,0) = (0,-1,0)
        let p = c.point_at(0.0);
        assert!((p - DVec3::new(0.0, -1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn offset_point_at_nonzero() {
        let c = make_offset_line();
        let p = c.point_at(5.0);
        assert!((p - DVec3::new(5.0, -1.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// GeomEval_SineWaveCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod sinewave_tests {
    use super::*;

    fn make_sine() -> Curve3 {
        Curve3::SineWave(SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        })
    }

    #[test]
    fn sine_at_zero() {
        let c = make_sine();
        let p = c.point_at(0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn sine_at_quarter_period() {
        let c = make_sine();
        let t = std::f64::consts::PI / 2.0;
        let p = c.point_at(t);
        // sin(pi/2) = 1, so y = 2.0
        assert!((p - DVec3::new(t, 2.0, 0.0)).length() < TOL);
    }

    #[test]
    fn sine_at_pi() {
        let c = make_sine();
        let p = c.point_at(std::f64::consts::PI);
        // sin(pi) = 0
        assert!((p - DVec3::new(std::f64::consts::PI, 0.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// GeomEval_CircularHelixCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod circular_helix_tests {
    use super::*;

    fn make_helix() -> Curve3 {
        Curve3::CircularHelix(CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 5.0,
            pitch: 2.0,
        })
    }

    #[test]
    fn helix_point_at_zero() {
        let h = make_helix();
        let p = h.point_at(0.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn helix_point_at_half_turn() {
        let h = make_helix();
        let p = h.point_at(std::f64::consts::PI);
        // After pi rad, point at (-R, 0, pitch/2)
        assert!((p - DVec3::new(-5.0, 0.0, 1.0)).length() < TOL);
    }

    #[test]
    fn helix_point_at_full_turn() {
        let h = make_helix();
        let p = h.point_at(std::f64::consts::TAU);
        // Full turn: back to (R,0,0) plus pitch
        assert!((p - DVec3::new(5.0, 0.0, 2.0)).length() < TOL);
    }
}

// =============================================================================
// GeomEval_EllipsoidSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod ellipsoid_surface_tests {
    use super::*;

    fn make_ellipsoid() -> Surface3 {
        Surface3::Ellipsoid(EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        })
    }

    #[test]
    fn ellipsoid_at_north_pole() {
        let s = make_ellipsoid();
        // v=0 → north pole: (0, 0, radius_z)
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::new(0.0, 0.0, 1.0)).length() < TOL);
    }

    #[test]
    fn ellipsoid_on_equator() {
        let s = make_ellipsoid();
        // u=0, v=pi/2 → equator along +X: (radius_x, 0, 0)
        let p = s.point_at(0.0, std::f64::consts::FRAC_PI_2);
        assert!((p - DVec3::new(3.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn ellipsoid_default_domain() {
        let s = make_ellipsoid();
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - std::f64::consts::TAU).abs() < TOL);
        assert!((v0 - 0.0).abs() < TOL);
        assert!((v1 - std::f64::consts::PI).abs() < TOL);
    }
}

// =============================================================================
// GeomEval_CircularHelicoidSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod helicoid_surface_tests {
    use super::*;

    fn make_helicoid() -> Surface3 {
        Surface3::Helicoid(HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 4.0,
        })
    }

    #[test]
    fn helicoid_at_u0v0() {
        let s = make_helicoid();
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn helicoid_default_domain() {
        let s = make_helicoid();
        let [u0, u1, v0, v1] = s.default_domain();
        let pi = std::f64::consts::PI;
        assert!((u0 - (-2.0 * pi)).abs() < TOL);
        assert!((u1 - (2.0 * pi)).abs() < TOL);
        assert!((v0 - (-10.0)).abs() < TOL);
        assert!((v1 - 10.0).abs() < TOL);
    }
}

// =============================================================================
// GeomAPI_ExtremaCurveCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod api_extrema_tests {
    use rcad_kernel::extrema_curve_curve;
    use super::*;

    #[test]
    fn parallel_lines_distance() {
        let l1 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let l2 = Curve3::Line(Line3 {
            origin: DVec3::new(0.0, 10.0, 0.0),
            direction: DVec3::X,
        });
        let ext = extrema_curve_curve(&l1, &l2, 32);
        assert!(ext.min_distance() - 10.0 < TOL, "distance={}", ext.min_distance());
    }

    #[test]
    fn intersecting_lines() {
        let l1 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let l2 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::Y });
        let ext = extrema_curve_curve(&l1, &l2, 32);
        assert!(ext.min_distance() < TOL, "skew distance={}", ext.min_distance());
    }
}

// =============================================================================
// GeomAPI_Interpolate_Test.cxx
// =============================================================================

#[cfg(test)]
mod api_interpolate_tests {
    use super::*;
    use rcad_kernel::fit::interpolate_points;

    #[test]
    fn interpolate_three_points() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(5.0, 0.0, 0.0),
            DVec3::new(10.0, 0.0, 0.0),
        ];
        let bs = interpolate_points(&pts).expect("interpolate");
        // At t=0 → first point
        let p0 = bs.point_at(0.0);
        assert!((p0 - pts[0]).length() < 1e-4, "first point mismatch");
        // At t=1 → last point
        let p1 = bs.point_at(1.0);
        assert!((p1 - pts[2]).length() < 1e-4, "last point mismatch");
    }
}
