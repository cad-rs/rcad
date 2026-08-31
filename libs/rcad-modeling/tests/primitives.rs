//! Integration tests for `rcad-modeling` primitive builders.
//!
//! Trimmed to the API that exists in rcad-modeling today (make_box_brep);
//! the sweep/loft/fillet tests were written against a planned API surface
//! that never landed and blocked `cargo test --workspace` from compiling.

use glam::DVec3;
use rcad_kernel::topods;
use rcad_kernel::topods::TShape;
use rcad_modeling::make_box_brep;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn face_count(brep: &topods::BRep) -> usize {
    brep.tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Face(_)))
        .count()
}

#[test]
fn box_face_count() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
    assert_eq!(face_count(&brep), 6, "box must have 6 faces");
}

#[test]
fn non_finite_dimension_returns_error() {
    let result = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, f64::NAN, 1.0, 1.0);
    assert!(
        result.is_err(),
        "NaN dimension should fail: {:?}",
        result.ok()
    );
}

