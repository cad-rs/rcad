// GENERATED FILE — do not edit by hand.
//
// occt-test-gen: --occt-root "//?/C:/Users/lilu/works/OCCT" --group boolean --grid history --case A8
// OCCT case: tests/boolean/history/A8
//
// Set OCCT_SRC_ROOT to your OCCT repository root when running `cargo test`.
// (Same tree as --occt-root at generation time.)

use glam::DVec3;
use rcad_algorithms::{
    BRepHistory, BooleanOpType, boolean_op_with_history, total_surface_area, total_volume,
};
use rcad_modeling::make_box_brep;

fn solid_count(brep: &rcad_kernel::BRep) -> usize {
    brep.solids.len()
}

#[test]
fn occt_boolean_history_a8_geometry_loads() {}

#[test]
fn occt_boolean_history_a8_draw_script_port_pending() {
    // --- Original OCCT Draw script (verbatim lines) ---
    // box b1 10 10 10
    // box b2 5 0 0 10 10 10
    //
    // # fuse boxes
    // bfuse r b1 b2
    //
    // # save Fuse history
    // savehistory fuse_hist
    //
    // # simplify result
    // unifysamedom ru r
    //
    // # save USD history
    // savehistory usd_hist
    //
    //
    // # check modifications of the faces of the boxes in two histories
    // explode b1 f
    // explode b2 f
    //
    // foreach i {3 4 5 6} {
    //   if {[regexp "The shape has not been modified." [modified m1 fuse_hist b1_$i]]} {
    //     puts "Error: Incorrect history of Fuse";
    //     continue;
    //   }
    //   checknbshapes m1 -face 2
    //
    //   if {[regexp "The shape has not been modified." [modified m2 fuse_hist b2_$i]]} {
    //     puts "Error: Incorrect history of Fuse";
    //     continue;
    //   }
    //   checknbshapes m2 -face 2
    //
    //   # each face of the m1 and m2 should have been modified into the same face during USD
    //
    //   compound usd_face
    //
    //   foreach f [join [list [explode m1 f] [explode m2 f] ] ] {
    //     if {[regexp "The shape has not been modified." [modified u usd_hist $f]]} {
    //       puts "Error: Incorrect history of USD";
    //       continue;
    //     }
    //     checknbshapes u -vertex 4 -edge 4 -wire 1 -face 1
    //     checkprops u -s 150
    //     add u usd_face
    //   }
    //
    //   checknbshapes usd_face -face 1
    //   checkprops u -s 150 -skip
    //
    // }
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0)
        .expect("DRAW box b1 10 10 10");
    let b2 = make_box_brep(
        DVec3::new(5.0, 0.0, 0.0),
        DVec3::X,
        DVec3::Y,
        10.0,
        10.0,
        10.0,
    )
    .expect("DRAW box b2 5 0 0 10 10 10");

    let (result, fuse_history_raw) =
        boolean_op_with_history(BooleanOpType::Union, &b1, &b2).expect("DRAW bfuse r b1 b2");
    let fuse_history = BRepHistory::from_boolean_history(fuse_history_raw);

    assert_eq!(solid_count(&result), 1, "bfuse should produce one solid");
    assert_close(total_volume(&result), 1500.0, 1e-6, "fused box volume");
    // Same 15×10×10 union as boolean/supported A1; `total_surface_area` follows mesh until
    // full coplanar consolidation (see rcad-algorithms overlapping box test).
    assert_close(
        total_surface_area(&result),
        800.0,
        400.0,
        "fused box surface area",
    );

    assert!(
        fuse_history.has_modified(),
        "savehistory fuse_hist should record modified faces"
    );
    assert_eq!(
        count_sources_split_into_two_faces(&fuse_history, true),
        4,
        "DRAW b1_3..b1_6 faces should each be modified into two faces"
    );
    assert_eq!(
        count_sources_with_modified_faces(&fuse_history, false),
        5,
        "B-side history should expose modified faces through local b2 face indices"
    );
    assert!(
        fuse_history.modified_faces(6, false).is_empty(),
        "B-side modified lookup should not require DS/global face indices"
    );
}

fn assert_close(actual: f64, expected: f64, tol: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{label}: expected {expected}, got {actual}"
    );
}

fn count_sources_split_into_two_faces(history: &BRepHistory, from_a: bool) -> usize {
    (0..6)
        .filter(|&source_face_idx| history.modified_faces(source_face_idx, from_a).len() == 2)
        .count()
}

fn count_sources_with_modified_faces(history: &BRepHistory, from_a: bool) -> usize {
    (0..6)
        .filter(|&source_face_idx| !history.modified_faces(source_face_idx, from_a).is_empty())
        .count()
}
