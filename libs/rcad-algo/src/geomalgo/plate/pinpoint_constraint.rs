//! OCCT Plate_PinpointConstraint (TKGeomAlgo/Plate) — 1:1 port.
//!
//! Plate_PinpointConstraint.hxx L28-53, .cxx L76-93, .lxx L110-128.
//! gp_XY -> DVec2, gp_XYZ -> DVec3 (architecture mapping).

use glam::{DVec2, DVec3};

/// OCCT Plate_PinpointConstraint — imposes a 3D value at a 2D point with
/// derivation orders (idu, idv).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PinpointConstraint {
    value: DVec3,
    pnt2d: DVec2,
    idu: i32,
    idv: i32,
}

impl Default for PinpointConstraint {
    /// OCCT Plate_PinpointConstraint() (.cxx L76-82).
    fn default() -> Self {
        PinpointConstraint {
            pnt2d: DVec2::ZERO,
            value: DVec3::ZERO,
            idu: 0,
            idv: 0,
        }
    }
}

impl PinpointConstraint {
    /// OCCT Plate_PinpointConstraint(point2d, ImposedValue, iu = 0, iv = 0)
    /// (.cxx L84-93).
    pub fn new(point2d: DVec2, imposed_value: DVec3, iu: i32, iv: i32) -> Self {
        PinpointConstraint {
            pnt2d: point2d,
            value: imposed_value,
            idu: iu,
            idv: iv,
        }
    }

    /// OCCT Pnt2d() (.lxx L110-113).
    pub fn pnt2d(&self) -> DVec2 {
        self.pnt2d
    }

    /// OCCT Idu() (.lxx L115-118).
    pub fn idu(&self) -> i32 {
        self.idu
    }

    /// OCCT Idv() (.lxx L120-123).
    pub fn idv(&self) -> i32 {
        self.idv
    }

    /// OCCT Value() (.lxx L125-128).
    pub fn value(&self) -> DVec3 {
        self.value
    }
}
