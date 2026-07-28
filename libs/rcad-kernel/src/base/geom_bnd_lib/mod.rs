//! Add geometry to bounding boxes (GeomBndLib).
//!
//! OCCT TKGeomBase GeomBndLib package.
//! Adds curves/surfaces to existing Bnd_Box with proper gap handling.

use crate::geom::{Curve3, Surface3};
use crate::math::bnd::BndBox;

/// Add a curve to a bounding box.
///
/// OCCT: `GeomBndLib::Add(Curve, Tol, Box)`.
pub fn add_curve_to_box(curve: &Curve3, box_: &mut BndBox, tol: f64) {
    if let Some([min, max]) = crate::base::bnd_lib::curve_bounding_box(curve) {
        box_.add_point(min);
        box_.add_point(max);
        box_.set_gap(tol);
    }
}

/// Add a surface to a bounding box.
///
/// OCCT: `GeomBndLib::Add(Surface, Tol, Box)`.
pub fn add_surface_to_box(surface: &Surface3, vertices: &[crate::topo::topology::Vertex], box_: &mut BndBox, tol: f64) {
    if let Some([min, max]) = crate::base::bnd_lib::surface_bounding_box(surface, vertices) {
        box_.add_point(min);
        box_.add_point(max);
        box_.set_gap(tol);
    }
}
