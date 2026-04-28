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
