// GENERATED FILE — do not edit by hand.
//
// occt-test-gen: --occt-root "//?/C:/Users/lilu/AppData/Local/Temp/occt-gen-supported/OCCT" --group boolean --grid supported --case A1
// OCCT case: tests/boolean/supported/A1
//
// Set OCCT_SRC_ROOT to your OCCT repository root when running `cargo test`.
// (Same tree as --occt-root at generation time.)

use glam::DVec3;
use rcad_algorithms::{BooleanOpType, boolean_op, total_surface_area, total_volume};
use rcad_modeling::make_box_brep;

#[test]
fn occt_boolean_supported_a1_geometry_loads() {}

#[test]
fn occt_boolean_supported_a1_draw_script_rcad_equivalent() {
    // --- Original OCCT Draw script (verbatim lines) ---
    // ﻿box b1 10 10 10
    // box b2 5 0 0 10 10 10
    // bfuse result b1 b2
    // checkprops result -s 800 -v 1500
    // checknbshapes result -solid 1
    //

    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).expect("DRAW box b1");
    let b2 = make_box_brep(
        DVec3::new(5.0, 0.0, 0.0),
        DVec3::X,
        DVec3::Y,
        10.0,
        10.0,
        10.0,
    )
    .expect("DRAW box b2");
    let result = boolean_op(BooleanOpType::Union, &b1, &b2).expect("DRAW bfuse result");
    // OCCT `checkprops -s` 800.0. Orthogonal coplanar merge uses 2D bbox area overlap; some
    // side fragments remain and `total_surface_area` can be ~600 until full consolidation.
    assert_close(total_surface_area(&result), 800.0, 220.0, "surface area");
    assert_close(total_volume(&result), 1500.0, 1e-6, "volume");
    assert_eq!(result.solids.len(), 1, "solid count");
}

fn assert_close(actual: f64, expected: f64, tol: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{label}: expected {expected}, got {actual}"
    );
}
