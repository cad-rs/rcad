//! OCCT GeomAdaptor — curve adaptor pattern with degenerated curve handling.
//!
//! OCCT source: src/ModelingData/TKG3d/GeomAdaptor/
//!
//! Provides a unified interface for curve evaluation with:
//! - Load(curve, first, last) with parameter validation
//! - Degenerated curve support (parameter range zero or within confusion)
//! - TransformedCurve with affine transform composition

use glam::{DVec3, DAffine3};
use rcad_kernel::geom::*;

const TOL: f64 = 1e-10;

/// OCCT Precision::Confusion() equivalent.
const CONFUSION: f64 = 1e-7;

/// Curve adaptor wrapping a Curve3 with parameter range and degenerated-curve support.
///
/// OCCT: GeomAdaptor_Curve
#[derive(Debug, Clone)]
pub struct CurveAdaptor {
    curve: Option<Curve3>,
    first: f64,
    last: f64,
    is_degenerated: bool,
    curve_type: Curve3Type,
}

/// Simplified curve type enum matching GeomAbs_CurveType.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Curve3Type {
    Line, Circle, Ellipse, Hyperbola, Parabola, BezierCurve, BSplineCurve,
    OffsetCurve, TrimmedCurve, OtherCurve,
}

impl From<&Curve3> for Curve3Type {
    fn from(c: &Curve3) -> Self {
        match c {
            Curve3::Line(_) => Curve3Type::Line,
            Curve3::Circle(_) => Curve3Type::Circle,
            Curve3::Ellipse(_) => Curve3Type::Ellipse,
            Curve3::Hyperbola(_) => Curve3Type::Hyperbola,
            Curve3::Parabola(_) => Curve3Type::Parabola,
            Curve3::Bezier(_) => Curve3Type::BezierCurve,
            Curve3::BSpline(_) => Curve3Type::BSplineCurve,
            Curve3::Offset(_) => Curve3Type::OffsetCurve,
            _ => Curve3Type::OtherCurve,
        }
    }
}

/// Error for invalid adaptor operations.
#[derive(Debug, Clone, PartialEq)]
pub enum AdaptorError {
    NullCurve,
    InvalidRange { first: f64, last: f64, confusion: f64 },
}

impl CurveAdaptor {
    pub fn new() -> Self {
        Self { curve: None, first: 0.0, last: 0.0, is_degenerated: false, curve_type: Curve3Type::OtherCurve }
    }

    /// OCCT: GeomAdaptor_Curve::Load(curve, UFirst, ULast)
    /// Validates that parameters are within confusion tolerance.
    pub fn load(&mut self, curve: Curve3, u_first: f64, u_last: f64) -> Result<(), AdaptorError> {
        // OCCT: if curve is null, throw Standard_NullObject
        // (Rust uses Result instead)

        // OCCT: check degeneracy — |last - first| <= Precision::Confusion()
        let diff = u_last - u_first;
        if diff.abs() <= CONFUSION {
            self.is_degenerated = true;
        } else if diff < -CONFUSION {
            // First > Last beyond confusion → error
            return Err(AdaptorError::InvalidRange {
                first: u_first, last: u_last, confusion: CONFUSION,
            });
        }

        let typ = Curve3Type::from(&curve);
        self.curve = Some(curve);
        self.first = u_first;
        self.last = u_last;
        self.curve_type = typ;
        Ok(())
    }

    /// OCCT: GeomAdaptor_Curve::Load(curve) — uses curve's own domain.
    pub fn load_default(&mut self, curve: Curve3) {
        let domain = curve.default_domain();
        let u_first = domain[0];
        let u_last = domain[1];
        let typ = Curve3Type::from(&curve);
        self.curve = Some(curve);
        self.first = u_first;
        self.last = u_last;
        self.curve_type = typ;
        self.is_degenerated = (u_last - u_first).abs() <= CONFUSION;
    }

    pub fn is_null(&self) -> bool { self.curve.is_none() }
    pub fn first_parameter(&self) -> f64 { self.first }
    pub fn last_parameter(&self) -> f64 { self.last }
    pub fn get_type(&self) -> Curve3Type { self.curve_type }
    pub fn is_periodic(&self) -> bool {
        matches!(&self.curve, Some(Curve3::Circle(_)) | Some(Curve3::Ellipse(_)))
    }
    pub fn is_closed(&self) -> bool {
        self.is_periodic()
    }

    /// OCCT: Value(u) — evaluate at parameter, clamped to valid range.
    pub fn value(&self, u: f64) -> Option<DVec3> {
        let curve = self.curve.as_ref()?;
        let u_clamped = if self.is_degenerated { self.first } else { u.clamp(self.first, self.last) };
        Some(curve.point_at(u_clamped))
    }

    /// OCCT: D1(u) — evaluate point + first derivative.
    pub fn d1(&self, u: f64) -> Option<(DVec3, DVec3)> {
        let curve = self.curve.as_ref()?;
        let u_clamped = if self.is_degenerated { self.first } else { u.clamp(self.first, self.last) };
        let pt = curve.point_at(u_clamped);
        let eps = 1e-7;
        let d1 = (curve.point_at(u_clamped + eps) - curve.point_at(u_clamped - eps)) / (2.0 * eps);
        Some((pt, d1))
    }
}

/// Transformed curve wrapping a Curve3 with an affine transform.
///
/// OCCT: GeomAdaptor_TransformedCurve
#[derive(Debug, Clone)]
pub struct TransformedCurveAdaptor {
    curve: Option<Curve3>,
    trsf: DAffine3,
    first: f64,
    last: f64,
    curve_type: Curve3Type,
}

impl TransformedCurveAdaptor {
    pub fn new() -> Self {
        Self { curve: None, trsf: DAffine3::IDENTITY, first: 0.0, last: 0.0, curve_type: Curve3Type::OtherCurve }
    }

    /// Construct with curve + transform (OCCT: GeomAdaptor_TransformedCurve(curve, trsf)).
    pub fn new_with_trsf(curve: Curve3, trsf: DAffine3) -> Self {
        let domain = curve.default_domain();
        let typ = Curve3Type::from(&curve);
        Self { curve: Some(curve), trsf, first: domain[0], last: domain[1], curve_type: typ }
    }

    /// Construct with curve + bounds + transform.
    pub fn new_with_bounds(curve: Curve3, first: f64, last: f64, trsf: DAffine3) -> Self {
        let typ = Curve3Type::from(&curve);
        Self { curve: Some(curve), trsf, first, last, curve_type: typ }
    }

    pub fn is_3d_curve(&self) -> bool { self.curve.is_some() }
    pub fn get_type(&self) -> Curve3Type { self.curve_type }
    pub fn trsf(&self) -> DAffine3 { self.trsf }
    pub fn first_parameter(&self) -> f64 { self.first }
    pub fn last_parameter(&self) -> f64 { self.last }

    /// Evaluate at parameter (applies transform).
    pub fn value(&self, u: f64) -> Option<DVec3> {
        let curve = self.curve.as_ref()?;
        let pt = curve.point_at(u);
        Some(self.trsf.transform_point3(pt))
    }

    /// Evaluate point + D1 (vectors transformed by rotation, not translation).
    pub fn d1(&self, u: f64) -> Option<(DVec3, DVec3)> {
        let curve = self.curve.as_ref()?;
        let pt = curve.point_at(u);
        let eps = 1e-7;
        let d1 = (curve.point_at(u + eps) - curve.point_at(u - eps)) / (2.0 * eps);
        // Transform point with full transform, D1 vector with rotation+scale only
        let trsf_pt = self.trsf.transform_point3(pt);
        let trsf_d1 = self.trsf.transform_vector3(d1);
        Some((trsf_pt, trsf_d1))
    }
}

// =============================================================================
// Tests — GeomAdaptor GTests
// =============================================================================

#[cfg(test)]
mod adaptor_curve_tests {
    use super::*;

    fn make_line() -> Curve3 {
        Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::new(1.0, 1.0, 0.0).normalize() })
    }

    fn make_circle() -> Curve3 {
        Curve3::Circle(Circle3::new(DVec3::new(5.0, 5.0, 0.0), DVec3::Z, 3.0))
    }

    #[test]
    fn load_valid_parameters() {
        let mut adaptor = CurveAdaptor::new();
        let result = adaptor.load(make_line(), 0.0, 10.0);
        assert!(result.is_ok());
        assert!((adaptor.first_parameter() - 0.0).abs() < TOL);
        assert!((adaptor.last_parameter() - 10.0).abs() < TOL);
    }

    #[test]
    fn load_equal_parameters_degenerated() {
        let mut adaptor = CurveAdaptor::new();
        assert!(adaptor.load(make_line(), 5.0, 5.0).is_ok());
        // Degenerated curve evaluates to a single point
        let pt = adaptor.value(5.0).unwrap();
        let expected = DVec3::ZERO + DVec3::new(1.0, 1.0, 0.0).normalize() * 5.0;
        assert!((pt - expected).length() < TOL);
    }

    #[test]
    fn load_parameters_within_confusion() {
        let mut adaptor = CurveAdaptor::new();
        assert!(adaptor.load(make_line(), 5.0, 5.0 + CONFUSION * 0.5).is_ok());
        assert!((adaptor.first_parameter() - 5.0).abs() < TOL);
    }

    #[test]
    fn load_first_greater_than_last_throws() {
        let mut adaptor = CurveAdaptor::new();
        let result = adaptor.load(make_line(), 10.0, 5.0);
        assert!(result.is_err());
        match result {
            Err(AdaptorError::InvalidRange { .. }) => {}
            _ => panic!("Expected InvalidRange error"),
        }
    }

    #[test]
    fn load_default_uses_curve_domain() {
        let mut adaptor = CurveAdaptor::new();
        adaptor.load_default(make_circle());
        assert!((adaptor.first_parameter() - 0.0).abs() < TOL);
        assert!((adaptor.last_parameter() - 2.0 * std::f64::consts::PI).abs() < TOL);
        assert!(adaptor.is_periodic());
    }

    #[test]
    fn degenerated_circle_at_single_point() {
        let mut adaptor = CurveAdaptor::new();
        assert!(adaptor.load(make_circle(), std::f64::consts::PI, std::f64::consts::PI).is_ok());
        let pt = adaptor.value(std::f64::consts::PI).unwrap();
        // Circle3::new: normal=Z, center=(5,5,0), radius=3
        // x_dir=Y, y_dir=-X. P(PI) = center + 3*(cos(PI)*Y + sin(PI)*(-X))
        // = (5,5,0) + 3*(-1)*Y = (5,5,0) + (0,-3,0) = (5,2,0)
        assert!((pt - DVec3::new(5.0, 2.0, 0.0)).length() < TOL);
    }
}

#[cfg(test)]
mod transformed_curve_tests {
    use super::*;

    fn translation(v: DVec3) -> DAffine3 {
        DAffine3::from_translation(v)
    }

    fn rotation_z(angle: f64) -> DAffine3 {
        let rot = glam::DQuat::from_rotation_z(angle);
        DAffine3::from_quat(rot)
    }

    #[test]
    fn default_constructor() {
        let tc = TransformedCurveAdaptor::new();
        assert!(tc.is_3d_curve() == false);
    }

    #[test]
    fn construct_with_curve_and_trsf() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let tc = TransformedCurveAdaptor::new_with_trsf(line, translation(DVec3::new(0.0, 0.0, 5.0)));
        assert!(tc.is_3d_curve());
        assert_eq!(tc.get_type(), Curve3Type::Line);
    }

    #[test]
    fn construct_with_bounds() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let tc = TransformedCurveAdaptor::new_with_bounds(line, 0.0, 10.0, translation(DVec3::new(0.0, 0.0, 5.0)));
        assert!((tc.first_parameter() - 0.0).abs() < TOL);
        assert!((tc.last_parameter() - 10.0).abs() < TOL);
    }

    #[test]
    fn line_identity_transform() {
        let line = Curve3::Line(Line3 { origin: DVec3::new(1.0, 2.0, 3.0), direction: DVec3::X });
        let tc = TransformedCurveAdaptor::new_with_trsf(line.clone(), DAffine3::IDENTITY);
        let pt = tc.value(5.0).unwrap();
        let expected = line.point_at(5.0);
        assert!((pt - expected).length() < TOL);
    }

    #[test]
    fn circle_identity_transform() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let curve = Curve3::Circle(circle);
        let tc = TransformedCurveAdaptor::new_with_trsf(curve.clone(), DAffine3::IDENTITY);
        for t in [0.0, 1.57, 3.14159, 4.71, 6.283] {
            let pt = tc.value(t).unwrap();
            let expected = curve.point_at(t);
            assert!((pt - expected).length() < TOL);
        }
    }

    #[test]
    fn line_translation() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let tc = TransformedCurveAdaptor::new_with_trsf(line.clone(), translation(DVec3::new(0.0, 0.0, 5.0)));
        let pt = tc.value(3.0).unwrap();
        assert!((pt - DVec3::new(3.0, 0.0, 5.0)).length() < TOL);
    }

    #[test]
    fn line_rotation() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let tc = TransformedCurveAdaptor::new_with_trsf(line.clone(), rotation_z(std::f64::consts::PI / 2.0));
        let pt = tc.value(3.0).unwrap();
        assert!((pt - DVec3::new(0.0, 3.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_translation() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let curve = Curve3::Circle(circle);
        let tc = TransformedCurveAdaptor::new_with_trsf(curve.clone(), translation(DVec3::new(10.0, 20.0, 30.0)));
        for t in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            let pt = tc.value(t).unwrap();
            let local = curve.point_at(t);
            assert!((pt - (local + DVec3::new(10.0, 20.0, 30.0))).length() < TOL);
        }
    }

    #[test]
    fn d1_translation() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let curve = Curve3::Circle(circle);
        let tc = TransformedCurveAdaptor::new_with_trsf(curve.clone(), translation(DVec3::new(10.0, 20.0, 30.0)));
        let (pt, d1) = tc.d1(0.0).unwrap();
        let local_pt = curve.point_at(0.0);
        assert!((pt - (local_pt + DVec3::new(10.0, 20.0, 30.0))).length() < TOL);
        // D1: translation doesn't change vectors
        let (_, local_d1) = CurveAdaptor { curve: Some(curve.clone()), first: 0.0, last: 6.283, is_degenerated: false, curve_type: Curve3Type::Circle }.d1(0.0).unwrap();
        assert!((d1 - local_d1).length() < TOL);
    }

    #[test]
    fn d1_rotation() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let tc = TransformedCurveAdaptor::new_with_trsf(line.clone(), rotation_z(std::f64::consts::PI / 2.0));
        let (_pt, d1) = tc.d1(3.0).unwrap();
        // Line D1 = (1,0,0), rotated 90deg → (0,1,0)
        assert!((d1 - DVec3::new(0.0, 1.0, 0.0)).length() < TOL);
    }
}
