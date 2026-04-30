//! Properties on **modeling-built** BReps: same invariants as kernel primitives, different construction path.
//!
//! Complements in-crate `PrimitiveSolid` tests by going through `rcad_modeling::make_box_brep`.

use rcad_algorithms::tolerance::*;
use glam::DVec3;
use rcad_algorithms::{boolean_op, total_surface_area, total_volume, BooleanOpType};
use rcad_kernel::face_surface_area;
use rcad_modeling::{make_box_brep, make_cone_brep, make_cylinder_brep};

#[test]
fn modeling_box_2x3x4_total_surface_area_and_volume() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");
    assert!((total_surface_area(&b) - 52.0).abs() < TOLERANCE_ADAPTIVE_MAX, "SA = 2(6+12+8)");
    assert!((total_volume(&b) - 24.0).abs() < TOLERANCE_ADAPTIVE_MAX, "V = 2*3*4");
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
    assert!((sum - total).abs() < TOLERANCE_ADAPTIVE_MAX, "per-face sum {sum} vs total {total}");
}

/// OCCT `bcommon_simple/B1`: 1³ ∩ 0.5×1×0.5 corner box → six faces, `checkprops -s 2.5`.
#[test]
fn occt_style_boolean_bcommon_simple_b1_surface_area() {
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let b1 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b1");
    let b2 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).expect("b2");
    let r = boolean_op(BooleanOpType::Intersection, &b1, &b2).expect("bcommon");
    let nf = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    assert_eq!(nf, 6, "expected six faces, got {nf}, area={a}");
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 2.5_f64);
    assert!(
        (a - 2.5).abs() <= tol,
        "surface area: expected 2.5, got {a} ({nf} faces)"
    );
}

/// OCCT `bcommon_simple/B3`: 1³ ∩ 0.5×1.5×1 (b2 at y −0.5); `checkprops -s 4`.
#[test]
fn occt_style_boolean_bcommon_simple_b3_surface_area() {
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let b1 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b1");
    let b2 = make_box_brep(DVec3::new(0.0, -0.5, 0.0), DVec3::X, DVec3::Y, 0.5, 1.5, 1.0).expect("b2");
    let r = boolean_op(BooleanOpType::Intersection, &b1, &b2).expect("bcommon");
    let nf = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    assert_eq!(nf, 6, "expected six faces, got {nf}, area={a}");
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 4.0_f64);
    assert!(
        (a - 4.0).abs() <= tol,
        "surface area: expected 4.0, got {a} ({nf} faces)"
    );
}

/// OCCT `bcommon_simple/C1`: 1³ ∩ 1.5×0.5×0.5; `checkprops -s 2.5`.
#[test]
fn occt_style_boolean_bcommon_simple_c1_surface_area() {
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let b1 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b1");
    let b2 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.5, 0.5, 0.5).expect("b2");
    let r = boolean_op(BooleanOpType::Intersection, &b1, &b2).expect("bcommon");
    let nf = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    assert_eq!(nf, 6, "expected six faces, got {nf}, area={a}");
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 2.5_f64);
    assert!(
        (a - 2.5).abs() <= tol,
        "surface area: expected 2.5, got {a} ({nf} faces)"
    );
}

/// OCCT `bcommon_simple/C3`: 1³ ∩ (0.5×0.5×1 at x=0.25); `checkprops -s 2.5`.
#[test]
fn occt_style_boolean_bcommon_simple_c3_surface_area() {
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let b1 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b1");
    let b2 = make_box_brep(DVec3::new(0.25, 0.0, 0.0), DVec3::X, DVec3::Y, 0.5, 0.5, 1.0).expect("b2");
    let r = boolean_op(BooleanOpType::Intersection, &b1, &b2).expect("bcommon");
    let nf = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    assert_eq!(nf, 6, "expected six faces, got {nf}, area={a}");
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 2.5_f64);
    assert!(
        (a - 2.5).abs() <= tol,
        "surface area: expected 2.5, got {a} ({nf} faces)"
    );
}

/// OCCT `bcommon_simple/A6`: two identical unit boxes, `bcommon` → full cube; `checkprops -s 6`.
#[test]
fn occt_style_boolean_bcommon_identical_unit_boxes_surface_area() {
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let b1 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b1");
    let b2 = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b2");
    let r = boolean_op(BooleanOpType::Intersection, &b1, &b2).expect("bcommon");
    let nf = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    assert_eq!(nf, 6, "expected six faces, got {nf}, area={a}");
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 6.0_f64);
    assert!(
        (a - 6.0).abs() <= tol,
        "surface area: expected 6.0, got {a} ({nf} faces)"
    );
}

/// OCCT `bcommon_simple/E1`: two upright prisms from z=0, half-width second box; `checkprops -s 4`.
#[test]
fn occt_style_boolean_bcommon_simple_e1_surface_area() {
    use glam::DVec3;
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let ba = make_box_brep(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0).normalize(),
        DVec3::new(0.0, 1.0, 0.0).normalize(),
        1.0,
        1.0,
        1.0,
    )
    .expect("ba");
    let bb = make_box_brep(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0).normalize(),
        DVec3::new(0.0, 1.0, 0.0).normalize(),
        1.0,
        0.5,
        1.0,
    )
    .expect("bb");
    let r = boolean_op(BooleanOpType::Intersection, &ba, &bb).expect("bcommon");
    let nf = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    assert_eq!(nf, 6, "expected six faces, got {nf}, area={a}");
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 4.0_f64);
    assert!(
        (a - 4.0).abs() <= tol,
        "surface area: expected 4.0, got {a} ({nf} faces)"
    );
}

/// OCCT `bcommon_simple/E3`: 2×2×2 box ∩ unit cube (common); `checkprops -s 6`.
#[test]
fn occt_style_boolean_bcommon_simple_e3_surface_area() {
    use glam::DVec3;
    use rcad_algorithms::{boolean_op, BooleanOpType};
    let ba = make_box_brep(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0).normalize(),
        DVec3::new(0.0, 1.0, 0.0).normalize(),
        2.0,
        2.0,
        2.0,
    )
    .expect("ba");
    let bb = make_box_brep(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0).normalize(),
        DVec3::new(0.0, 1.0, 0.0).normalize(),
        1.0,
        1.0,
        1.0,
    )
    .expect("bb");
    let r = boolean_op(BooleanOpType::Intersection, &ba, &bb).expect("bcommon");
    let nf = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    assert_eq!(nf, 6, "expected six faces, got {nf}, area={a}");
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 6.0_f64);
    assert!(
        (a - 6.0).abs() <= tol,
        "surface area: expected 6.0, got {a} ({nf} faces)"
    );
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
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 4.41421_f64);
    assert!(
        (area - 4.41421).abs() <= tol,
        "surface area: expected ~4.41421, got {area}"
    );
}

/// OCCT `bopcommon_simple/ZP7`: `pcone` ∩ prism cylinder on exploded base face (`checkprops -s 919.56`).
#[test]
fn occt_style_bopcommon_simple_zp7_cone_cylinder_intersection_surface_area() {
    let pc = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 10.0, 20.0).expect("cone");
    let pcy = make_cylinder_brep(
        DVec3::new(0.0, 0.0, -5.0),
        DVec3::Z,
        DVec3::X,
        10.0,
        10.0,
    )
    .expect("cylinder");
    let r = boolean_op(BooleanOpType::Intersection, &pc, &pcy).expect("intersection");
    let nf: usize = r
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let a = total_surface_area(&r);
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 919.56_f64);
    assert!(
        (a - 919.56).abs() <= tol,
        "surface area: expected ~919.56, got {a} ({nf} faces, {} solids)",
        r.solids.len()
    );
}

/// OCCT `boptuc_simple/ZP3`: `boptuc` is `Cut(G2,G1)` → cylinder − cone (`checkprops -s 1390.8`).
///
/// Expected to hit the coaxial `cylinder \\ cone` loft-shell shortcut inside `boolean_op`.
#[test]
fn occt_style_boptuc_simple_zp3_cylinder_minus_cone_surface_area() {
    let pc = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 10.0, 20.0).expect("cone");
    let pcy = make_cylinder_brep(
        DVec3::new(0.0, 0.0, -5.0),
        DVec3::Z,
        DVec3::X,
        10.0,
        10.0,
    )
    .expect("cylinder");
    let r = boolean_op(BooleanOpType::Difference, &pcy, &pc).expect("cylinder minus cone");
    let a = total_surface_area(&r);
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 1390.8_f64);
    assert!(
        (a - 1390.8).abs() <= tol,
        "surface area: expected ~1390.8, got {a}"
    );
}

/// OCCT `bopcut_simple/ZP8`: sharp cone − prism cylinder on exploded base (`checkprops -s 254.16`).
#[test]
fn occt_style_bopcut_simple_zp8_cone_minus_cylinder_surface_area() {
    let pc = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 10.0, 20.0).expect("cone");
    let pcy = make_cylinder_brep(
        DVec3::new(0.0, 0.0, -5.0),
        DVec3::Z,
        DVec3::X,
        10.0,
        10.0,
    )
    .expect("cylinder");
    let r = boolean_op(BooleanOpType::Difference, &pc, &pcy).expect("difference");
    let a = total_surface_area(&r);
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 254.16_f64);
    assert!(
        (a - 254.16).abs() <= tol,
        "surface area: expected ~254.16, got {a}"
    );
}

/// OCCT `bfuse_simple/E5`: offset box ∪ extruded unit square prism (`checkprops -s 170`).
#[test]
fn occt_style_bfuse_simple_e5_box_union_offset_prism_surface_area() {
    let ba = make_box_brep(DVec3::new(3.0, 3.0, 0.0), DVec3::X, DVec3::Y, 5.0, 7.0, 4.0)
        .expect("ba");
    let bb = make_box_brep(DVec3::new(3.0, 2.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("bb");
    let r = boolean_op(BooleanOpType::Union, &ba, &bb).expect("union");
    let a = total_surface_area(&r);
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 170.0_f64);
    assert!(
        (a - 170.0).abs() <= tol,
        "surface area: expected ~170 (OCCT checkprops -s), got {a}"
    );
}
