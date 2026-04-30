//! FEA / mechanical multi-scale workflow presets and fuzzy reporting.

use rcad_algorithms::tolerance::*;
use glam::DVec3;
use rcad_algorithms::{
    boolean_op_robust, boolean_op_with_options, resolved_boolean_fuzzy_tol_for_ds, tolerance,
    BooleanOpType, BooleanOptions, BooleanRobustOptions, ExtremeGeometryRetryPolicy,
};
use rcad_modeling::make_box_brep;

#[test]
fn resolved_fuzzy_clamps_to_confusion_floor() {
    assert_eq!(
        resolved_boolean_fuzzy_tol_for_ds(0.0),
        tolerance::TOLERANCE_ABS
    );
    assert_eq!(
        resolved_boolean_fuzzy_tol_for_ds(-1.0),
        tolerance::TOLERANCE_ABS
    );
    assert_eq!(
        resolved_boolean_fuzzy_tol_for_ds(TOLERANCE_MESH_LEGACY),
        TOLERANCE_MESH_LEGACY
    );
}

#[test]
fn bare_options_zero_fuzzy_reports_effective_floor() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("a");
    let b = make_box_brep(DVec3::new(0.2, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b");
    let opts = BooleanOptions {
        fuzzy_tol: 0.0,
        ..Default::default()
    };
    let (_, report) = boolean_op_with_options(BooleanOpType::Intersection, &a, &b, opts).expect("op");
    assert_eq!(report.configured_fuzzy_tol, 0.0);
    assert_eq!(
        report.effective_fuzzy_tol,
        tolerance::TOLERANCE_ABS
    );
}

#[test]
fn propagate_geom_tolerances_flag_runs_when_enabled() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("a");
    let b = make_box_brep(DVec3::new(0.2, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b");
    let opts = BooleanOptions {
        run_propagate_geom_tolerances: true,
        fuzzy_tol: 0.0,
        ..Default::default()
    };
    let (out, report) =
        boolean_op_with_options(BooleanOpType::Intersection, &a, &b, opts).expect("op");
    assert!(report.propagated_geom_tolerances);
    assert!(
        out.geom.vertex_tolerance.len() >= out.vertices.len(),
        "propagate_tolerances sizes vertex_tolerance"
    );
}

#[test]
fn fea_preset_enables_glue_healing_make_connected_and_scaled_fuzzy() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
    let b = make_box_brep(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");
    let ro = BooleanRobustOptions::for_fea(&a, &b);
    assert!(ro.base.use_glue);
    assert!(ro.base.run_healing);
    assert!(ro.base.run_make_connected);
    assert!(ro.base.run_propagate_geom_tolerances);
    assert!(ro.base.fuzzy_tol >= tolerance::TOLERANCE_ABS);
    assert!(ro.base.glue_tolerance >= tolerance::TOLERANCE_ABS);
    assert_eq!(ro.fuzzy_retry_ladder.len(), 3);
    assert!(ro.fuzzy_retry_ladder[0] > tolerance::TOLERANCE_ABS);

    let (_out, report) = boolean_op_robust(
        BooleanOpType::Intersection,
        &a,
        &b,
        ro,
    )
    .expect("robust fea preset");
    assert!(report.effective_fuzzy_tol >= tolerance::TOLERANCE_ABS);
    assert!(
        report.propagated_geom_tolerances,
        "FEA preset runs propagate_tolerances after pipeline"
    );
}

#[test]
fn mechanical_multiscale_preset_is_geometry_aware_with_extended_ladder() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("a");
    let b = make_box_brep(
        DVec3::new(500.0, 0.0, 0.0),
        DVec3::X,
        DVec3::Y,
        1.0,
        1.0,
        1.0,
    )
    .expect("b");

    let ro = BooleanRobustOptions::for_mechanical_multiscale(&a, &b);
    assert!(
        matches!(
            ro.extreme_geometry.policy,
            ExtremeGeometryRetryPolicy::GeometryAware
        ),
        "expected GeometryAware policy"
    );
    assert!(ro.fuzzy_retry_ladder.len() >= 5);
    assert!(ro.base.use_glue);
    assert!(ro.base.run_propagate_geom_tolerances);
    assert!(ro.base.fuzzy_tol >= tolerance::TOLERANCE_ABS);
}
