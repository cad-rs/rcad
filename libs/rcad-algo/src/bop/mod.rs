//! TKBO ?Boolean Operation algorithms (OCCT TKBO toolkit).

pub mod algo;
pub mod ds;
pub mod tools;
pub mod brep_algo_api;
pub mod int_tools;

use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Surface3};

/// Minimal 3D point classifier.
pub fn classify_point(_point: DVec3, _face_indices: &[usize], _ds: &ds::DS) -> ds::Classification {
    ds::Classification::Out
}

/// 3D curve projection wrapper — routes through topalgo.
pub fn closest_point_on_curve(curve: &Curve3, query: DVec3) -> (f64, DVec3) {
    crate::topalgo::adaptor::project_point_on_curve(curve, query)
}

/// 3D surface projection wrapper — routes through topalgo.
pub fn closest_point_on_surface(surface: &Surface3, point: DVec3) -> (DVec2, DVec3) {
    crate::topalgo::adaptor::project_point_on_surface(surface, point)
}

/// Curve bounding box (delegates to rcad-kernel BndLib, no tolerance param).
pub fn curve_bounds(curve: &Curve3) -> (DVec3, DVec3) {
    match rcad_kernel::base::bnd_lib::curve_bounding_box(curve) {
        Some([min, max]) => (min, max),
        None => (DVec3::ZERO, DVec3::ZERO),
    }
}

/// Curve bounding box with range (delegates to rcad-kernel BndLib).
pub fn curve_bounds_with_range(curve: &Curve3, _first: f64, _last: f64) -> (DVec3, DVec3) {
    match rcad_kernel::base::bnd_lib::curve_bounding_box(curve) {
        Some([min, max]) => (min, max),
        None => (DVec3::ZERO, DVec3::ZERO),
    }
}
