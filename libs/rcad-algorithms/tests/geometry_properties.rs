//! Properties on **modeling-built** BReps: same invariants as kernel primitives, different construction path.
//!
//! Complements in-crate `PrimitiveSolid` tests by going through `rcad_modeling::make_box_brep`.

use glam::DVec3;
use rcad_algorithms::{total_surface_area, total_volume};
use rcad_kernel::face_surface_area;
use rcad_modeling::make_box_brep;

#[test]
fn modeling_box_2x3x4_total_surface_area_and_volume() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");
    assert!((total_surface_area(&b) - 52.0).abs() < 1e-3, "SA = 2(6+12+8)");
    assert!((total_volume(&b) - 24.0).abs() < 1e-3, "V = 2*3*4");
}

#[test]
fn modeling_box_face_surface_areas_sum_to_total() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");
    let total = total_surface_area(&b);
    let mut sum = 0.0;
    let mut i = 0usize;
    for solid in &b.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                sum += face_surface_area(&b, face, i);
                i += 1;
            }
        }
    }
    assert!((sum - total).abs() < 1e-3, "per-face sum {sum} vs total {total}");
}

/// OCCT `bcommon_simple/C8`: two boxes, one rotated; `checkprops -s 4.41421`.
#[test]
fn occt_style_boolean_bcommon_simple_c8_surface_area() {
    use glam::{DAffine3, DVec3};
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b1");
    let mut b2 = make_box_brep(
        DVec3::ZERO,
        DVec3::X,
        DVec3::Y,
        std::f64::consts::SQRT_2,
        std::f64::consts::FRAC_1_SQRT_2,
        1.0,
    )
    .expect("b2");
    let pivot = DVec3::new(0.0, 0.0, 0.0);
    let axis = DVec3::new(0.0, 0.0, 1.0).normalize_or(DVec3::Z);
    let rot = DAffine3::from_axis_angle(axis, (45.0_f64).to_radians());
    let xf = DAffine3::from_translation(pivot) * rot * DAffine3::from_translation(-pivot);
    b2.apply_transform(xf);
    let r = boolean_op(BooleanOpType::Intersection, &b1, &b2).expect("intersection");
    let area = total_surface_area(&r);
    let tol = (5e-3_f64).max(0.0625 * 4.41421_f64);
    assert!(
        (area - 4.41421).abs() <= tol,
        "surface area: expected ~4.41421, got {area}"
    );
}
