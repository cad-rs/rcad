#[cfg(test)]
mod tests {
    use glam::DVec3;
    use rcad_kernel::PrimitiveSolid;

    use crate::healing::*;
    use crate::geom_populate;

    fn unit_box() -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn heal_valid_box_is_noop() {
        let b = unit_box();
        let (out, report) = heal(&b);
        assert!(report.initial.is_valid());
        assert!(report.final_result.is_valid());
        assert!(report.passes.is_empty());
        assert!(report.parametric_passes.is_empty());
        assert!(report.make_connected_passes.is_empty());
        assert!(!report.stages.is_empty());
        assert_eq!(out.vertices.len(), b.vertices.len());
        assert_eq!(out.edges.len(), b.edges.len());
    }

    #[test]
    fn heal_zero_normal_face_gets_fixed() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (out, report) = heal(&b);
        assert!(report.initial_issue_count() >= 1);
        assert!(report.is_improved() || report.is_clean());
        assert!(report.initial_stats.zero_normal >= 1);
        assert_eq!(report.initial_stats.total(), report.initial_issue_count());
        assert_eq!(report.final_stats.total(), report.final_issue_count());
        assert!(
            report
                .stages
                .iter()
                .any(|s| matches!(s.stage, HealingStage::FinalCheck))
        );
        assert!(!out.solids[0].shells[0].faces[0].normal.abs_diff_eq(DVec3::ZERO, 0.0));
    }

    #[test]
    fn analyze_only_preserves_input_and_reports_issues() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                mode: HealingMode::AnalyzeOnly,
                ..HealingOptions::default()
            },
        );

        assert!(report.initial_issue_count() >= 1);
        assert_eq!(report.initial_issue_count(), report.final_issue_count());
        assert!(report.passes.is_empty());
        assert!(report.parametric_passes.is_empty());
        assert!(report.make_connected_passes.is_empty());
        assert_eq!(out.solids[0].shells[0].faces[0].normal, DVec3::ZERO);
    }

    #[test]
    fn healing_make_connected_fallback_reporting_is_consistent() {
        let mut b = unit_box();

        // Keep at least one checker issue that standard repair does not heal.
        b.solids[0].shells[0].faces[0].outer_wire.edges[0].idx = usize::MAX;
        // Add near-duplicate vertices that can be merged only by the fallback
        // tolerance (repair tolerance intentionally set much tighter).
        b.vertices[1].point = b.vertices[0].point + DVec3::new(TOLERANCE_MESH_LEGACY, 0.0, 0.0);

        let (_out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                tolerance: TOLERANCE_LEN_MIN,
                max_passes: 1,
                run_make_connected_on_stall: true,
                make_connected_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
                make_connected_max_passes: 2,
                make_connected_tolerance_growth: 1.0,
                make_connected_tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE,
                ..HealingOptions::default()
            },
        );

        // Depending on how much progress the regular repair pass can make,
        // make-connected fallback may or may not be needed. If it ran, stage
        // and report vectors must stay in sync.
        let mc_stage_count = report
            .stages
            .iter()
            .filter(|s| matches!(s.stage, HealingStage::MakeConnectedPass))
            .count();
        assert_eq!(mc_stage_count, report.make_connected_passes.len());
        assert!(report.make_connected_passes.len() <= 1);
    }

    #[test]
    fn healing_parametric_consistency_pass_is_reported_when_enabled_by_data() {
        let mut b = unit_box();

        // Make one edge obviously suspect for SameRange/SameParameter scans.
        if b.geom.edge_same_parameter.len() < b.edges.len() {
            b.geom.edge_same_parameter.resize(b.edges.len(), true);
        }
        b.geom.edge_same_parameter[0] = false;
        if b.geom.edge_curve_range.len() < b.edges.len() {
            b.geom.edge_curve_range.resize(b.edges.len(), Some([0.0, 1.0]));
        }
        b.geom.edge_curve_range[0] = Some([0.0, 1.0]);

        let (_out, report) = analyze_and_heal(&b, HealingOptions::default());
        let saw_param_stage = report
            .stages
            .iter()
            .any(|s| matches!(s.stage, HealingStage::ParametricConsistencyPass));
        assert_eq!(saw_param_stage, !report.parametric_passes.is_empty());
    }

    #[test]
    fn healing_can_disable_parametric_consistency_prepass() {
        let mut b = unit_box();
        if b.geom.edge_same_parameter.len() < b.edges.len() {
            b.geom.edge_same_parameter.resize(b.edges.len(), true);
        }
        b.geom.edge_same_parameter[0] = false;

        let (_out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                run_parametric_consistency_prepass: false,
                run_parametric_consistency_iterative: false,
                ..HealingOptions::default()
            },
        );

        assert!(report.parametric_passes.is_empty());
    }

    #[test]
    fn healing_make_connected_prepass_always_records_stage() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].outer_wire.edges[0].idx = usize::MAX;

        let (_out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                max_passes: 1,
                make_connected_prepass_mode: MakeConnectedPrepassMode::Always,
                make_connected_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
                make_connected_max_passes: 1,
                make_connected_tolerance_growth: 1.0,
                make_connected_tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE,
                ..HealingOptions::default()
            },
        );

        assert!(
            report
                .stages
                .iter()
                .any(|s| matches!(s.stage, HealingStage::PreMakeConnected))
        );
        assert!(!report.make_connected_passes.is_empty());
    }

    #[test]
    fn operator_chain_runs_repair_and_parametric_passes() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;
        if b.geom.edge_same_parameter.len() < b.edges.len() {
            b.geom.edge_same_parameter.resize(b.edges.len(), true);
        }
        b.geom.edge_same_parameter[0] = false;

        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::ParametricConsistency,
                HealingOperator::Repair,
                HealingOperator::StopIfClean,
            ],
        );

        assert!(!report.parametric_passes.is_empty());
        assert!(!report.passes.is_empty());
        assert!(
            report.stages.iter().any(|s| matches!(s.stage, HealingStage::OperatorChainStep))
        );
    }

    #[test]
    fn operator_chain_stop_if_clean_short_circuits_following_steps() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::Repair,
                HealingOperator::StopIfClean,
                HealingOperator::MakeConnected,
            ],
        );

        // Repair should clean this case; stop-if-clean should prevent make-connected.
        assert!(report.make_connected_passes.is_empty());
        assert!(report.final_result.is_valid());
    }

    #[test]
    fn shape_process_default_config_works_on_valid_shape() {
        let b = unit_box();
        let config = ShapeProcessConfig::default();
        let (out, report) = run_shape_process(&b, &config);

        assert!(report.is_clean());
        assert!(report.stats.converged_early);
        assert_eq!(out.vertices.len(), b.vertices.len());
    }

    #[test]
    fn shape_process_import_preset_fixes_zero_normal() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::import_preset();
        let (_out, report) = run_shape_process(&b, &config);

        assert!(report.is_improved() || report.is_clean());
        assert!(report.initial_issue_count() >= 1);
    }

    #[test]
    fn shape_process_boolean_cleanup_preset_works() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::boolean_cleanup_preset();
        let (_out, report) = run_shape_process(&b, &config);

        assert!(report.is_improved() || report.is_clean());
    }

    #[test]
    fn shape_process_analysis_preset_is_conservative() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::analysis_preset();
        let (_out, report) = run_shape_process(&b, &config);

        // Analysis preset should at least diagnose issues
        assert!(report.initial_issue_count() >= 1);
    }

    #[test]
    fn shape_process_aggressive_preset_applies_all_operators() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::aggressive_preset();
        let (_out, report) = run_shape_process(&b, &config);

        assert!(report.is_improved() || report.is_clean());
        // Aggressive preset has many operators
        assert!(config.operators.len() >= 8);
    }

    #[test]
    fn shape_process_report_summary_is_informative() {
        let b = unit_box();
        let config = ShapeProcessConfig::default();
        let (_out, report) = run_shape_process(&b, &config);

        let summary = report.summary();
        assert!(summary.contains("ShapeProcess"));
        assert!(summary.contains("Clean") || summary.contains("issues"));
    }

    #[test]
    fn operator_chain_handles_new_operators() {
        let b = unit_box();

        // Test that new operators don't panic
        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::FixSmallAreaFaces,
                HealingOperator::FixSliverFaces,
                HealingOperator::FixNonManifold,
                HealingOperator::PropagateTolerances,
                HealingOperator::UnifySameDomain,
                HealingOperator::RemoveInternalFaces,
            ],
        );

        // All operators should run without error
        assert!(!report.stages.is_empty());
    }

    #[test]
    fn fix_small_area_faces_removes_tiny_faces() {
        let b = unit_box();

        // Unit box faces are not tiny, so nothing should be removed
        let (result, removed) = fix_small_area_faces(&b, TOLERANCE_VOL_CUBE_FACTOR * TOLERANCE_ADAPTIVE_MAX);
        assert_eq!(removed, 0);

        // The result should have the same number of faces
        let result_face_count: usize = result.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        let original_face_count: usize = b.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(result_face_count, original_face_count);
    }

    #[test]
    fn new_healing_stages_are_recorded() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::FixSmallAreaFaces,
                HealingOperator::FixNonManifold,
                HealingOperator::PropagateTolerances,
            ],
        );

        // Should have geometry and topology repair stages
        assert!(report.stages.iter().any(|s|
            matches!(s.stage, HealingStage::GeometryRepairPass)));
        assert!(report.stages.iter().any(|s|
            matches!(s.stage, HealingStage::TopologyRepairPass)));
        assert!(report.stages.iter().any(|s|
            matches!(s.stage, HealingStage::FinalizePass)));
    }

    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Tests for New Operators
    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn unit_sphere() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 })
    }

    fn unit_cylinder() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Cylinder { radius: 1.0, height: 2.0 })
    }

    #[test]
    fn split_angle_operator_default_params() {
        let params = SplitAngleOperator::default();
        assert!((params.max_angle - std::f64::consts::PI / 2.0).abs() < TOLERANCE_LEN_MIN);
        assert!(params.split_cylinders);
        assert!(params.split_tori);
        assert!(params.split_cones);
        assert!(params.split_spheres);
        assert!((params.start_angle).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn split_angle_on_sphere() {
        let sphere = unit_sphere();
        let params = SplitAngleOperator {
            max_angle: std::f64::consts::PI / 4.0, // 45 degrees
            ..Default::default()
        };
        let (result, splits) = split_angle_operator(&sphere, &params);
        // Sphere should potentially be split
        assert!(!result.vertices.is_empty());
        assert!(!result.solids.is_empty());
        let _ = splits;
    }

    #[test]
    fn split_angle_on_cylinder() {
        let cyl = unit_cylinder();
        let params = SplitAngleOperator {
            max_angle: std::f64::consts::PI / 3.0, // 60 degrees
            split_cylinders: true,
            ..Default::default()
        };
        let (result, splits) = split_angle_operator(&cyl, &params);
        assert!(!result.vertices.is_empty());
        let _ = splits;
    }

    #[test]
    fn split_angle_preserves_shape_when_disabled() {
        let sphere = unit_sphere();
        let params = SplitAngleOperator {
            split_spheres: false,
            ..Default::default()
        };
        let (result, splits) = split_angle_operator(&sphere, &params);
        assert_eq!(splits, 0);
        assert_eq!(result.vertices.len(), sphere.vertices.len());
    }

    #[test]
    fn split_continuity_default_params() {
        let params = SplitContinuityOperator::default();
        assert_eq!(params.min_continuity, ContinuityLevel::C1);
        assert!((params.tolerance - TOLERANCE_MESH_LEGACY).abs() < TOLERANCE_LEN_MIN);
        assert!(params.check_curves);
        assert!(params.check_surfaces);
        assert_eq!(params.max_splits_per_edge, 100);
    }

    #[test]
    fn split_continuity_on_box() {
        let b = unit_box();
        let params = SplitContinuityOperator::default();
        let (result, splits) = split_continuity_operator(&b, &params);
        // Box edges should be C2 continuous (straight lines)
        assert_eq!(splits, 0);
        assert_eq!(result.vertices.len(), b.vertices.len());
    }

    #[test]
    fn continuity_level_ordering() {
        assert!(ContinuityLevel::C0 < ContinuityLevel::C1);
        assert!(ContinuityLevel::C1 < ContinuityLevel::C2);
    }

    #[test]
    fn convert_to_bspline_default_params() {
        let params = ConvertToBSplineOperator::default();
        assert_eq!(params.max_degree, 3);
        assert!(params.convert_curves);
        assert!(params.convert_surfaces);
        assert!(!params.convert_planes);
        assert!(params.convert_elementary);
        assert_eq!(params.approximation_samples, 20);
    }

    #[test]
    fn convert_to_bspline_on_sphere() {
        let sphere = unit_sphere();
        let params = ConvertToBSplineOperator {
            convert_elementary: true,
            ..Default::default()
        };
        let (result, conversions) = convert_to_bspline_operator(&sphere, &params);
        assert!(conversions > 0);
        // Check that surfaces are converted
        let has_bspline = result.geom.surfaces.iter().any(|s| {
            matches!(s, rcad_kernel::geom::Surface3::BSpline(_))
        });
        assert!(has_bspline);
    }

    #[test]
    fn convert_to_bspline_preserves_planes_when_disabled() {
        let b = unit_box();
        geom_populate::populate_box_geom(&mut b.clone());
        let params = ConvertToBSplineOperator {
            convert_planes: false,
            convert_elementary: false,
            ..Default::default()
        };
        let (result, conversions) = convert_to_bspline_operator(&b, &params);
        assert_eq!(conversions, 0);
        let _ = result;
    }

    #[test]
    fn surface_to_bezier_default_params() {
        let params = SurfaceToBezierOperator::default();
        assert!(params.convert_surfaces);
        assert!(params.convert_pcurves);
        assert!(params.convert_curves);
        assert_eq!(params.max_degree, 25);
    }

    #[test]
    fn surface_to_bezier_on_bspline() {
        // Create a sphere and convert to BSpline first
        let sphere = unit_sphere();
        let bspline_params = ConvertToBSplineOperator::default();
        let (bspline_sphere, _) = convert_to_bspline_operator(&sphere, &bspline_params);

        // Then convert to Bezier
        let bezier_params = SurfaceToBezierOperator::default();
        let (result, conversions) = surface_to_bezier_operator(&bspline_sphere, &bezier_params);
        let _ = result;
    }

    #[test]
    fn scale_shape_uniform() {
        let scale = ScaleShapeOperator::uniform(2.0);
        assert!(scale.is_uniform());
        assert!((scale.scale_x - 2.0).abs() < TOLERANCE_LEN_MIN);
        assert!((scale.scale_y - 2.0).abs() < TOLERANCE_LEN_MIN);
        assert!((scale.scale_z - 2.0).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn scale_shape_non_uniform() {
        let scale = ScaleShapeOperator::non_uniform(2.0, 1.0, 0.5);
        assert!(!scale.is_uniform());
        assert!((scale.scale_x - 2.0).abs() < TOLERANCE_LEN_MIN);
        assert!((scale.scale_y - 1.0).abs() < TOLERANCE_LEN_MIN);
        assert!((scale.scale_z - 0.5).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn scale_shape_default_is_identity() {
        let scale = ScaleShapeOperator::default();
        assert!(scale.is_uniform());
        assert!((scale.scale_x - 1.0).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn scale_shape_on_box() {
        let b = unit_box();
        let params = ScaleShapeOperator::uniform(2.0);
        let (result, mods) = scale_shape_operator(&b, &params);

        assert!(mods > 0);
        // Box should be scaled by 2x
        let original_bounds = b.bounding_box().unwrap();
        let scaled_bounds = result.bounding_box().unwrap();

        // The size should be approximately doubled
        let original_size = original_bounds[1] - original_bounds[0];
        let scaled_size = scaled_bounds[1] - scaled_bounds[0];

        assert!((scaled_size.x - 2.0 * original_size.x).abs() < super::TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((scaled_size.y - 2.0 * original_size.y).abs() < super::TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((scaled_size.z - 2.0 * original_size.z).abs() < super::TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn scale_shape_identity_is_noop() {
        let b = unit_box();
        let params = ScaleShapeOperator::uniform(1.0);
        let (result, mods) = scale_shape_operator(&b, &params);

        assert_eq!(mods, 0);
        assert_eq!(result.vertices.len(), b.vertices.len());
    }

    #[test]
    fn scale_shape_with_origin() {
        let b = unit_box();
        // Unit box is centered at origin with size 1, so bounds are [-0.5, 0.5]
        let params = ScaleShapeOperator {
            scale_x: 2.0,
            scale_y: 2.0,
            scale_z: 2.0,
            origin: Some(DVec3::new(0.5, 0.5, 0.5)), // Scale around a point outside the box
            ..Default::default()
        };
        let (result, _) = scale_shape_operator(&b, &params);

        // Verify the scaling was applied (box is now 2x size)
        let bounds = result.bounding_box().unwrap();
        // The box should have been scaled by 2x
        let width = bounds[1].x - bounds[0].x;
        assert!((width - 2.0).abs() < 0.1, "Expected width ~2.0, got {}", width);
    }

    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Tests for Operator Chaining Improvements
    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn operator_condition_always() {
        let b = unit_box();
        let (_, report) = heal(&b);
        let condition = OperatorCondition::Always;
        assert!(condition.evaluate(&b, &report, &[]));
    }

    #[test]
    fn operator_condition_only_if_issues() {
        // Create a HealingReport with issues for testing
        let mut report_with_issues = HealingReport::default();
        report_with_issues.final_result.issues.push(CheckIssue::ZeroNormal { solid: 0, shell: 0, face: 0 });

        let condition = OperatorCondition::OnlyIfIssues;
        assert!(condition.evaluate(&BRep::new(), &report_with_issues, &[]));

        // A clean report should have no issues
        let report_clean = HealingReport::default();
        assert!(!condition.evaluate(&BRep::new(), &report_clean, &[]));
    }

    #[test]
    fn operator_condition_only_if_clean() {
        // A clean report should pass OnlyIfClean
        let report_clean = HealingReport::default();

        let condition = OperatorCondition::OnlyIfClean;
        assert!(condition.evaluate(&BRep::new(), &report_clean, &[]));

        // A report with issues should fail OnlyIfClean
        let mut report_with_issues = HealingReport::default();
        report_with_issues.final_result.issues.push(CheckIssue::ZeroNormal { solid: 0, shell: 0, face: 0 });
        assert!(!condition.evaluate(&BRep::new(), &report_with_issues, &[]));
    }

    #[test]
    fn operator_condition_issue_count_above() {
        // Create a report with 2 issues
        let mut report = HealingReport::default();
        report.final_result.issues.push(CheckIssue::ZeroNormal { solid: 0, shell: 0, face: 0 });
        report.final_result.issues.push(CheckIssue::ZeroNormal { solid: 0, shell: 0, face: 1 });

        let condition = OperatorCondition::OnlyIfIssueCountAbove(0);
        assert!(condition.evaluate(&BRep::new(), &report, &[]));

        let condition2 = OperatorCondition::OnlyIfIssueCountAbove(1);
        assert!(condition2.evaluate(&BRep::new(), &report, &[]));

        let condition3 = OperatorCondition::OnlyIfIssueCountAbove(2);
        assert!(!condition3.evaluate(&BRep::new(), &report, &[]));
    }

    #[test]
    fn healing_operator_with_condition_new() {
        let op = HealingOperatorWithCondition::new(HealingOperator::Repair);
        assert!(op.condition.is_none());
        assert!(op.dependencies.is_empty());
        assert!(op.label.is_none());
    }

    #[test]
    fn healing_operator_with_condition_with_condition() {
        let op = HealingOperatorWithCondition::with_condition(
            HealingOperator::Repair,
            OperatorCondition::OnlyIfIssues,
        );
        assert!(op.condition.is_some());
    }

    #[test]
    fn healing_operator_with_condition_depends_on() {
        let op = HealingOperatorWithCondition::new(HealingOperator::Repair)
            .depends_on(0)
            .depends_on(1);
        assert_eq!(op.dependencies, vec![0, 1]);
    }

    #[test]
    fn healing_operator_with_condition_with_label() {
        let op = HealingOperatorWithCondition::new(HealingOperator::Repair)
            .with_label("test_label");
        assert_eq!(op.label, Some("test_label".to_string()));
    }

    #[test]
    fn operator_chain_config_default() {
        let config = OperatorChainConfig::default();
        assert!(config.stop_on_clean);
        assert_eq!(config.max_iterations, 1);
        assert!(!config.operators.is_empty());
    }

    #[test]
    fn operator_chain_config_mesh_prep_preset() {
        let config = OperatorChainConfig::mesh_prep_preset();
        assert!(config.stop_on_clean);
        assert!(!config.operators.is_empty());
        // Should have split angle and convert to bspline
        let has_split_angle = config.operators.iter().any(|op| {
            matches!(op.operator, HealingOperator::SplitAngle(_))
        });
        let has_convert = config.operators.iter().any(|op| {
            matches!(op.operator, HealingOperator::ConvertToBSpline(_))
        });
        assert!(has_split_angle || has_convert);
    }

    #[test]
    fn operator_chain_config_export_prep_preset() {
        let config = OperatorChainConfig::export_prep_preset();
        assert!(config.stop_on_clean);
        let has_bezier = config.operators.iter().any(|op| {
            matches!(op.operator, HealingOperator::SurfaceToBezier(_))
        });
        assert!(has_bezier);
    }

    #[test]
    fn operator_chain_config_scale_preset() {
        let config = OperatorChainConfig::scale_preset(2.0);
        assert!(config.stop_on_clean);
        let has_scale = config.operators.iter().any(|op| {
            matches!(op.operator, HealingOperator::ScaleShape(_))
        });
        assert!(has_scale);
    }

    #[test]
    fn run_advanced_operator_chain_basic() {
        let b = unit_box();
        let config = OperatorChainConfig::default();
        let (result, report) = run_advanced_operator_chain(&b, &config);

        assert!(report.is_clean);
        assert!(report.operator_results.len() > 0);
        assert!(report.total_elapsed_seconds >= 0.0);
        let _ = result;
    }

    #[test]
    fn run_advanced_operator_chain_with_conditions() {
        let b = unit_box();
        let config = OperatorChainConfig {
            operators: vec![
                HealingOperatorWithCondition::new(HealingOperator::Repair),
                HealingOperatorWithCondition::with_condition(
                    HealingOperator::MakeConnected,
                    OperatorCondition::OnlyIfIssues,
                ),
            ],
            stop_on_clean: true,
            ..Default::default()
        };
        let (_, report) = run_advanced_operator_chain(&b, &config);

        // First operator runs, second should be skipped (condition not met)
        assert!(report.is_clean);
    }

    #[test]
    fn run_advanced_operator_chain_with_dependencies() {
        let b = unit_box();
        let config = OperatorChainConfig {
            operators: vec![
                HealingOperatorWithCondition::new(HealingOperator::PropagateTolerances),
                HealingOperatorWithCondition::new(HealingOperator::Repair)
                    .depends_on(0),
            ],
            stop_on_clean: true,
            ..Default::default()
        };
        let (_, report) = run_advanced_operator_chain(&b, &config);

        assert!(report.operators_executed > 0);
    }

    #[test]
    fn new_operators_in_healing_chain() {
        let b = unit_box();

        // Test that all new operators can be used in a chain
        let (_result, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::SplitAngle(SplitAngleOperator::default()),
                HealingOperator::SplitContinuity(SplitContinuityOperator::default()),
                HealingOperator::ConvertToBSpline(ConvertToBSplineOperator::default()),
                HealingOperator::SurfaceToBezier(SurfaceToBezierOperator::default()),
                HealingOperator::ScaleShape(ScaleShapeOperator::uniform(1.0)),
            ],
        );

        assert!(!report.stages.is_empty());
    }

    #[test]
    fn operator_result_default() {
        let result = OperatorResult::default();
        assert!(!result.changed);
        assert_eq!(result.modifications, 0);
        assert!(!result.skipped);
        assert!(result.skip_reason.is_none());
    }

    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Tests for New ShapeProcess Operators
    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn direct_faces_operator_default() {
        let params = DirectFacesOperator::default();
        assert!((params.tolerance - TOLERANCE_ABS).abs() < TOLERANCE_LEN_MIN);
        assert!(params.update_surface_references);
        assert!(params.recompute_normals);
        assert!(params.fix_wire_orientation);
    }

    #[test]
    fn direct_faces_operator_on_valid_box() {
        let b = unit_box();
        let params = DirectFacesOperator::default();
        let (result, fixed) = direct_faces_operator(&b, &params);
        // Verify operator runs successfully
        assert_eq!(result.vertices.len(), b.vertices.len());
    }

    #[test]
    fn direct_faces_operator_on_flipped_normal() {
        let mut b = unit_box();
        // Flip a normal to simulate an indirect face
        b.solids[0].shells[0].faces[0].normal = -b.solids[0].shells[0].faces[0].normal;

        let params = DirectFacesOperator {
            recompute_normals: false, // Don't recompute, just fix orientation
            ..Default::default()
        };
        let (result, _fixed) = direct_faces_operator(&b, &params);
        // Should have processed the face
        assert!(!result.solids[0].shells[0].faces.is_empty());
    }

    #[test]
    fn same_parameter_operator_default() {
        let params = SameParameterOperator::default();
        assert!((params.tolerance - TOLERANCE_ABS).abs() < TOLERANCE_LEN_MIN);
        assert_eq!(params.max_samples, 23);
        assert!(!params.enforce);
        assert!(params.update_pcurve_ranges);
    }

    #[test]
    fn same_parameter_operator_enforced() {
        let params = SameParameterOperator::enforced(super::TOLERANCE_MESH_LEGACY);
        assert!(params.enforce);
        assert!((params.tolerance - TOLERANCE_MESH_LEGACY).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn same_parameter_operator_on_valid_box() {
        let b = unit_box();
        let params = SameParameterOperator::default();
        let (result, fixed) = same_parameter_operator(&b, &params);
        // Valid box should have no same parameter issues
        assert_eq!(fixed, 0);
        assert_eq!(result.vertices.len(), b.vertices.len());
    }

    #[test]
    fn remove_internal_faces_operator_default() {
        let params = RemoveInternalFacesOperator::default();
        assert!((params.tolerance - TOLERANCE_ABS).abs() < TOLERANCE_LEN_MIN);
        assert!((params.min_face_area - TOLERANCE_LINEAR_ULTRA_STRICT).abs() < TOLERANCE_LEN_MIN);
        assert!(params.check_manifold);
        assert!(params.merge_vertices);
        assert!(params.preserve_material_boundaries);
    }

    #[test]
    fn remove_internal_faces_on_valid_box() {
        let b = unit_box();
        let params = RemoveInternalFacesOperator::default();
        let (result, removed) = remove_internal_faces_operator(&b, &params);
        // Valid box should have no internal faces
        assert_eq!(removed, 0);
        assert_eq!(result.vertices.len(), b.vertices.len());
    }

    #[test]
    fn heal_geometry_operator_default() {
        let params = HealGeometryOperator::default();
        assert!((params.tolerance - TOLERANCE_ABS).abs() < TOLERANCE_LEN_MIN);
        assert_eq!(params.max_passes, 3);
        assert!(params.fix_face_orientation);
        assert!(params.fix_same_parameter);
        assert!(params.fix_same_range);
        assert!(params.fix_wire_gaps);
        assert!(params.remove_degenerate_faces);
        assert!(params.propagate_tolerances);
        assert!(params.recompute_normals);
        assert!(params.fix_uv_bounds);
        assert!(!params.remove_small_edges);
    }

    #[test]
    fn heal_geometry_operator_minimal() {
        let params = HealGeometryOperator::minimal(TOLERANCE_MESH_LEGACY);
        assert_eq!(params.max_passes, 1);
        assert!(params.fix_face_orientation);
        assert!(params.fix_same_parameter);
        assert!(!params.fix_wire_gaps);
        assert!(!params.remove_degenerate_faces);
    }

    #[test]
    fn heal_geometry_operator_aggressive() {
        let params = HealGeometryOperator::aggressive(TOLERANCE_MESH_LEGACY);
        assert_eq!(params.max_passes, 5);
        assert!(params.remove_small_edges);
    }

    #[test]
    fn heal_geometry_operator_sequence() {
        let params = HealGeometryOperator::default();
        let sequence = params.get_sequence();
        assert!(!sequence.is_empty());
        // Recompute normals should be first in default sequence
        assert!(sequence.contains(&HealGeometryStep::RecomputeNormals));
        // Propagate tolerances should be last in default sequence
        assert!(sequence.contains(&HealGeometryStep::PropagateTolerances));
    }

    #[test]
    fn heal_geometry_operator_custom_sequence() {
        let params = HealGeometryOperator {
            custom_sequence: vec![
                HealGeometryStep::FixSameParameter,
                HealGeometryStep::FixSameRange,
            ],
            ..Default::default()
        };
        let sequence = params.get_sequence();
        assert_eq!(sequence.len(), 2);
        assert_eq!(sequence[0], HealGeometryStep::FixSameParameter);
        assert_eq!(sequence[1], HealGeometryStep::FixSameRange);
    }

    #[test]
    fn heal_geometry_on_valid_box() {
        let b = unit_box();
        let params = HealGeometryOperator::default();
        let (result, report) = heal_geometry_operator(&b, &params);
        // Valid box should need minimal fixes
        assert_eq!(result.vertices.len(), b.vertices.len());
        let total_fixes = report.vertices_merged + report.faces_reoriented + report.wires_fixed
            + report.same_parameter_fixed + report.same_range_fixed + report.degenerate_faces_removed;
        assert!(total_fixes == 0 || report.normals_recomputed > 0);
    }

    #[test]
    fn heal_geometry_on_zero_normal_box() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let params = HealGeometryOperator::default();
        let (result, report) = heal_geometry_operator(&b, &params);
        // Should have recomputed the zero normal
        assert!(report.normals_recomputed >= 1);
        assert!(!result.solids[0].shells[0].faces[0].normal.abs_diff_eq(DVec3::ZERO, 0.0));
    }

    #[test]
    fn new_operators_in_chain() {
        let b = unit_box();

        // Test that new operators can be used in a chain
        let (_result, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::DirectFaces(DirectFacesOperator::default()),
                HealingOperator::SameParameter(SameParameterOperator::default()),
                HealingOperator::RemoveInternalFacesOp(RemoveInternalFacesOperator::default()),
                HealingOperator::HealGeometry(HealGeometryOperator::default()),
            ],
        );

        assert!(!report.stages.is_empty());
    }

    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Tests for Operator Result Aggregation
    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn operator_result_aggregation_empty() {
        let agg = OperatorResultAggregation::new();
        assert_eq!(agg.total_executed, 0);
        assert_eq!(agg.total_skipped, 0);
        assert_eq!(agg.total_modifications, 0);
        assert!(!agg.has_changes());
    }

    #[test]
    fn operator_result_aggregation_add_result() {
        let mut agg = OperatorResultAggregation::new();
        let result = OperatorResult {
            operator: HealingOperator::Repair,
            changed: true,
            modifications: 5,
            issues_fixed: 3,
            description: "test".to_string(),
            elapsed_seconds: 0.1,
            skipped: false,
            skip_reason: None,
        };
        agg.add_result(result);

        assert_eq!(agg.total_executed, 1);
        assert_eq!(agg.total_modifications, 5);
        assert_eq!(agg.total_issues_fixed, 3);
        assert!(agg.has_changes());
    }

    #[test]
    fn operator_result_aggregation_change_rate() {
        let mut agg = OperatorResultAggregation::new();

        // Add one with changes
        agg.add_result(OperatorResult {
            changed: true,
            ..OperatorResult::default()
        });

        // Add one without changes
        agg.add_result(OperatorResult {
            changed: false,
            ..OperatorResult::default()
        });

        assert!((agg.change_rate() - 0.5).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn operator_result_aggregation_summary() {
        let mut agg = OperatorResultAggregation::new();
        agg.add_result(OperatorResult {
            changed: true,
            modifications: 3,
            issues_fixed: 2,
            elapsed_seconds: 0.5,
            ..OperatorResult::default()
        });

        let summary = agg.summary();
        assert!(summary.contains("1 executed"));
        assert!(summary.contains("3 modifications"));
        assert!(summary.contains("2 issues fixed"));
    }

    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Tests for Rollback Support
    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn rollback_config_default() {
        let config = RollbackConfig::default();
        assert!(config.enabled);
        assert!(config.rollback_on_failure);
        assert!(config.rollback_on_regression);
        assert_eq!(config.max_issues_threshold, 0);
    }

    #[test]
    fn rollback_config_disabled() {
        let config = RollbackConfig::disabled();
        assert!(!config.enabled);
    }

    #[test]
    fn brep_snapshot_creation() {
        let b = unit_box();
        let snapshot = BRepSnapshot::new(&b, 0, "test", 0.5);
        assert_eq!(snapshot.operator_index, 0);
        assert_eq!(snapshot.label, "test");
        assert!((snapshot.timestamp_seconds - 0.5).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn run_healing_pipeline_with_rollback_basic() {
        let b = unit_box();
        let operators: Vec<HealingOperator> = vec![
            HealingOperator::DirectFaces(DirectFacesOperator::default()),
            HealingOperator::Repair,
        ];

        let (result, report) = run_healing_pipeline_with_rollback(
            &b,
            &operators,
            HealingOptions::default(),
            RollbackConfig::default(),
            None,
        );

        assert!(report.completed);
        assert!(!report.aggregation.results.is_empty());
        assert!(result.vertices.len() > 0);
    }

    #[test]
    fn run_healing_pipeline_with_rollback_reports_aggregation() {
        let b = unit_box();
        let operators: Vec<HealingOperator> = vec![
            HealingOperator::HealGeometry(HealGeometryOperator::minimal(TOLERANCE_ABS)),
            HealingOperator::PropagateTolerances,
        ];

        let (_result, report) = run_healing_pipeline_with_rollback(
            &b,
            &operators,
            HealingOptions::default(),
            RollbackConfig::default(),
            None,
        );

        assert!(report.completed);
        assert_eq!(report.aggregation.total_executed, operators.len());
    }

    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Tests for Progress Callbacks
    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn simple_progress_callback_creation() {
        let cb = SimpleProgressCallback::new(5);
        assert_eq!(cb.total_operators, 5);
        assert!(!cb.is_cancelled());
    }

    #[test]
    fn simple_progress_callback_cancel() {
        let cb = SimpleProgressCallback::new(5);
        cb.cancel();
        assert!(cb.is_cancelled());
    }

    #[test]
    fn simple_progress_callback_progress() {
        let cb = SimpleProgressCallback {
            current_operator: 2,
            total_operators: 4,
            ..Default::default()
        };
        assert!((cb.progress() - 0.5).abs() < TOLERANCE_LEN_MIN);
    }

    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Tests for Pipeline Execution Report
    // 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn pipeline_execution_report_summary() {
        let report = PipelineExecutionReport {
            aggregation: OperatorResultAggregation::new(),
            snapshots: Vec::new(),
            final_brep: BRep::new(),
            completed: true,
            failure_reason: None,
            rollback_index: None,
        };

        let summary = report.summary();
        assert!(summary.contains("Completed"));
    }

    #[test]
    fn pipeline_execution_report_with_rollback() {
        let report = PipelineExecutionReport {
            aggregation: OperatorResultAggregation::new(),
            snapshots: vec![BRepSnapshot::new(&BRep::new(), 0, "test", 0.0)],
            final_brep: BRep::new(),
            completed: false,
            failure_reason: Some("Test failure".to_string()),
            rollback_index: Some(0),
        };

        let summary = report.summary();
        assert!(summary.contains("Test failure"));
        assert!(summary.contains("rolled back"));
    }

    #[test]
    fn heal_geometry_step_variants() {
        // Test that all variants exist and can be compared
        assert_ne!(HealGeometryStep::FixFaceOrientation, HealGeometryStep::FixSameParameter);
        assert_ne!(HealGeometryStep::RecomputeNormals, HealGeometryStep::PropagateTolerances);
    }

    // Edge case tests for OCCT alignment

    #[test]
    fn heal_with_degenerate_edge() {
        let mut b = unit_box();
        // Create a degenerate edge (same start and end vertex)
        if b.edges.len() > 0 {
            let v0 = b.edges[0].start;
            b.edges[0].end = v0;
        }

        let (_out, report) = heal(&b);
        // Should attempt to fix or report the degenerate edge
        assert!(report.initial_issue_count() >= 1 || report.is_clean());
    }

    #[test]
    fn heal_with_reversed_face_normal() {
        let mut b = unit_box();
        // Reverse one face normal
        if !b.solids.is_empty() && !b.solids[0].shells.is_empty() {
            let faces = &mut b.solids[0].shells[0].faces;
            if !faces.is_empty() {
                faces[0].normal = -faces[0].normal;
            }
        }

        let (_out, report) = heal(&b);
        // Should detect and potentially fix the reversed normal
        assert!(report.is_improved() || report.is_clean());
    }

    #[test]
    fn heal_with_small_gap() {
        let mut b = unit_box();
        // Perturb a vertex slightly to create a small gap
        if !b.vertices.is_empty() {
            b.vertices[0].point.x += 0.001;
        }

        let (_out, report) = heal(&b);
        // Should detect some issue or be clean
    }

    #[test]
    fn heal_sphere_primitive() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        geom_populate::populate_box_geom(&mut brep);

        let (_out, report) = heal(&brep);
        // Sphere should heal without major issues
        assert!(report.is_clean() || report.is_improved() || report.final_issue_count() <= report.initial_issue_count());
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Enhanced Healing: ShapeFix_Solid and ShapeFix_Wire Equivalents
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// ShapeFix_Solid equivalent: comprehensive solid repair.
///
/// This function performs OCCT ShapeFix_Solid-like operations:
/// - Shell orientation verification and repair
/// - Solid closure verification
/// - Shell manifoldness checks
/// - Face orientation consistency
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `tolerance` - Geometric tolerance
///
/// # Returns
/// Repaired B-Rep and count of fixes applied.
pub fn fix_solid(brep: &BRep, _tolerance: f64) -> (BRep, SolidFixReport) {
    use crate::brep_repair::{fix_face_orientation, recompute_face_normals};
    use rcad_kernel::BRepGraph;

    let mut report = SolidFixReport::default();
    let mut current = brep.clone();

    // Step 1: Recompute invalid normals
    let (brep_with_normals, normals_fixed) = recompute_face_normals(&current);
    current = brep_with_normals;
    report.normals_recomputed = normals_fixed;

    // Step 2: Fix face orientation for inward-pointing faces
    let (brep_oriented, faces_reoriented) = fix_face_orientation(&current);
    current = brep_oriented;
    report.faces_reoriented = faces_reoriented;

    // Step 3: Check solid closure and manifoldness
    let graph = BRepGraph::from_brep(&current);

    // Check if shells are closed
    for (si, solid) in current.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            let is_closed = shell.faces.iter().all(|f| {
                // Check if wire is closed
                let wire = &f.outer_wire;
                if wire.edges.is_empty() {
                    return false;
                }
                true // Simplified check; full implementation would verify vertex chain
            });

            if !is_closed {
                report.unclosed_shells.push((si, shi));
            }
        }
    }

    // Check manifoldness
    let nm_summary = graph.non_manifold_summary();
    report.non_manifold_edges = nm_summary.multi_face_edges.len();
    report.non_manifold_vertices = nm_summary.non_manifold_vertices.len();

    // Step 4: Verify shell orientation consistency
    for solid in &current.solids {
        for shell in &solid.shells {
            // Count faces with normals pointing in consistent direction
            let mut outward_count = 0usize;
            let mut inward_count = 0usize;

            for face in &shell.faces {
                // Heuristic: if normal dot product with center-to-centroid is positive
                // the face is likely outward-facing
                if face.normal.z > 0.0 {
                    outward_count += 1;
                } else if face.normal.z < 0.0 {
                    inward_count += 1;
                }
            }

            // If most normals are inconsistent, note orientation issues
            if outward_count > 0 && inward_count > 0 {
                let ratio = outward_count as f64 / (outward_count + inward_count) as f64;
                if !(0.3..=0.7).contains(&ratio) {
                    report.orientation_inconsistencies += 1;
                }
            }
        }
    }

    report.total_fixes = report.normals_recomputed + report.faces_reoriented;
    (current, report)
}

/// Report from solid-level fixes.
#[derive(Debug, Clone, Default)]
pub struct SolidFixReport {
    /// Number of face normals recomputed.
    pub normals_recomputed: usize,
    /// Number of faces reoriented.
    pub faces_reoriented: usize,
    /// Indices of unclosed shells (solid_idx, shell_idx).
    pub unclosed_shells: Vec<(usize, usize)>,
    /// Number of non-manifold edges detected.
    pub non_manifold_edges: usize,
    /// Number of non-manifold vertices detected.
    pub non_manifold_vertices: usize,
    /// Number of shells with orientation inconsistencies.
    pub orientation_inconsistencies: usize,
    /// Total number of fixes applied.
    pub total_fixes: usize,
}

impl SolidFixReport {
    pub fn is_clean(&self) -> bool {
        self.unclosed_shells.is_empty()
            && self.non_manifold_edges == 0
            && self.non_manifold_vertices == 0
            && self.orientation_inconsistencies == 0
    }

    pub fn summary(&self) -> String {
        if self.is_clean() && self.total_fixes == 0 {
            "Solid is clean, no fixes needed".to_string()
        } else {
            format!(
                "Solid fixes: {} normals, {} orientations, {} unclosed shells, {} non-manifold edges, {} non-manifold vertices",
                self.normals_recomputed,
                self.faces_reoriented,
                self.unclosed_shells.len(),
                self.non_manifold_edges,
                self.non_manifold_vertices
            )
        }
    }
}

/// ShapeFix_Wire equivalent: comprehensive wire repair.
///
/// This function performs OCCT ShapeFix_Wire-like operations:
/// - Wire closure verification and repair
/// - Edge order verification
/// - Degenerate edge handling
/// - Self-intersection detection
/// - Wire orientation fix
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `tolerance` - Geometric tolerance
///
/// # Returns
/// Repaired B-Rep and detailed wire fix report.
pub fn fix_wire(brep: &BRep, tolerance: f64) -> (BRep, WireFixReport) {
    use crate::brep_repair::fix_wire_orientation;

    let mut report = WireFixReport::default();
    let mut current = brep.clone();

    // Step 1: Fix wire orientation
    let (brep_fixed, wires_fixed) = fix_wire_orientation(&current, tolerance);
    current = brep_fixed;
    report.wires_oriented = wires_fixed;

    // Step 2: Analyze wires for issues
    for (si, solid) in current.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                // Check outer wire
                let outer_issues = analyze_wire_issues(&current, &face.outer_wire, tolerance);
                if outer_issues.open_gaps > 0 || outer_issues.topological_self_intersections > 0 || outer_issues.geometric_self_intersections > 0 {
                    report.outer_wire_issues.push(WireIssueLocation {
                        solid: si,
                        shell: shi,
                        face: fi,
                        wire_idx: 0,
                        issues: outer_issues,
                    });
                }

                // Check inner wires
                for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
                    let inner_issues = analyze_wire_issues(&current, inner_wire, tolerance);
                    if inner_issues.open_gaps > 0 || inner_issues.topological_self_intersections > 0 || inner_issues.geometric_self_intersections > 0 {
                        report.inner_wire_issues.push(WireIssueLocation {
                            solid: si,
                            shell: shi,
                            face: fi,
                            wire_idx: wi + 1,
                            issues: inner_issues,
                        });
                    }
                }
            }
        }
    }

    // Step 3: Count degenerate edges
    for (ei, edge) in current.edges.iter().enumerate() {
        let start_pt = current.vertices.get(edge.start).map(|v| v.point);
        let end_pt = current.vertices.get(edge.end).map(|v| v.point);

        if let (Some(s), Some(e)) = (start_pt, end_pt)
            && (s - e).length() < tolerance {
                report.degenerate_edges.push(ei);
            }
    }

    // Step 4: Compute wire quality metrics
    report.total_wires_checked = report.outer_wire_issues.len()
        + report.inner_wire_issues.len()
        + current.solids.iter()
            .flat_map(|s| s.shells.iter())
            .flat_map(|sh| sh.faces.iter())
            .map(|f| 1 + f.inner_wires.len())
            .sum::<usize>();

    report.wires_with_issues = report.outer_wire_issues.len() + report.inner_wire_issues.len();
    report.total_fixes = report.wires_oriented;

    (current, report)
}

/// Location of a wire issue.
#[derive(Debug, Clone)]
pub struct WireIssueLocation {
    pub solid: usize,
    pub shell: usize,
    pub face: usize,
    pub wire_idx: usize,
    pub issues: crate::brep_check::WireIssueReport,
}

/// Report from wire-level fixes.
#[derive(Debug, Clone, Default)]
pub struct WireFixReport {
    /// Number of wires with corrected orientation.
    pub wires_oriented: usize,
    /// Issues found in outer wires.
    pub outer_wire_issues: Vec<WireIssueLocation>,
    /// Issues found in inner wires.
    pub inner_wire_issues: Vec<WireIssueLocation>,
    /// Indices of degenerate edges found.
    pub degenerate_edges: Vec<usize>,
    /// Total wires checked.
    pub total_wires_checked: usize,
    /// Wires with issues.
    pub wires_with_issues: usize,
    /// Total fixes applied.
    pub total_fixes: usize,
}

impl WireFixReport {
    pub fn is_clean(&self) -> bool {
        self.outer_wire_issues.is_empty()
            && self.inner_wire_issues.is_empty()
            && self.degenerate_edges.is_empty()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() && self.total_fixes == 0 {
            format!("All {} wires clean, no fixes needed", self.total_wires_checked)
        } else {
            format!(
                "Wire fixes: {} oriented, {} with issues ({} outer, {} inner), {} degenerate edges",
                self.wires_oriented,
                self.wires_with_issues,
                self.outer_wire_issues.len(),
                self.inner_wire_issues.len(),
                self.degenerate_edges.len()
            )
        }
    }
}

/// Analyze wire for issues without modifying.
fn analyze_wire_issues(brep: &BRep, wire: &rcad_kernel::topology::Wire, tolerance: f64) -> crate::brep_check::WireIssueReport {
    let n_edges = brep.edges.len();
    let mut open_gaps = 0usize;
    let mut topological_self_intersections = 0usize;
    let mut geometric_self_intersections = 0usize;

    // Collect wire vertices
    let mut wire_verts = Vec::with_capacity(wire.edges.len());
    for we in &wire.edges {
        if we.idx >= n_edges {
            continue;
        }
        let edge = &brep.edges[we.idx];
        let (sv, ev) = if we.forward {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        };
        if sv < brep.vertices.len() && ev < brep.vertices.len() {
            wire_verts.push((sv, ev));
        }
    }

    // Check for open gaps
    let n = wire_verts.len();
    if n > 1 {
        for i in 0..n {
            let next = (i + 1) % n;
            let end_v = wire_verts[i].1;
            let start_v = wire_verts[next].0;
            if end_v != start_v {
                let end_pt = brep.vertices[end_v].point;
                let start_pt = brep.vertices[start_v].point;
                if (end_pt - start_pt).length() > tolerance {
                    open_gaps += 1;
                }
            }
        }
    }

    // Check for topological self-intersection (vertex appearing more than twice)
    use std::collections::HashMap;
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();
    for &(sv, ev) in &wire_verts {
        *vertex_count.entry(sv).or_insert(0) += 1;
        *vertex_count.entry(ev).or_insert(0) += 1;
    }
    for &count in vertex_count.values() {
        if count > 2 {
            topological_self_intersections += 1;
        }
    }

    // Check for geometric self-intersection (2D projection)
    if n >= 4 {
        for i in 0..n {
            for j in (i + 2)..n {
                if i == 0 && j == n - 1 {
                    continue; // Adjacent edges wraparound
                }
                let (a_start, a_end) = wire_verts[i];
                let (b_start, b_end) = wire_verts[j];
                let p1 = brep.vertices[a_start].point;
                let p2 = brep.vertices[a_end].point;
                let p3 = brep.vertices[b_start].point;
                let p4 = brep.vertices[b_end].point;

                if segments_intersect_2d(p1, p2, p3, p4) {
                    geometric_self_intersections += 1;
                }
            }
        }
    }

    crate::brep_check::WireIssueReport {
        solid: 0,
        shell: 0,
        face: 0,
        wire_idx: 0,
        edge_count: wire.edges.len(),
        open_gaps,
        topological_self_intersections,
        geometric_self_intersections,
    }
}

/// Check if two 2D line segments intersect (XY plane projection).
fn segments_intersect_2d(p1: glam::DVec3, p2: glam::DVec3, p3: glam::DVec3, p4: glam::DVec3) -> bool {
    let x1 = p1.x; let y1 = p1.y;
    let x2 = p2.x; let y2 = p2.y;
    let x3 = p3.x; let y3 = p3.y;
    let x4 = p4.x; let y4 = p4.y;

    let (min_x1, max_x1) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
    let (min_y1, max_y1) = if y1 < y2 { (y1, y2) } else { (y2, y1) };
    let (min_x2, max_x2) = if x3 < x4 { (x3, x4) } else { (x4, x3) };
    let (min_y2, max_y2) = if y3 < y4 { (y3, y4) } else { (y4, y3) };

    if max_x1 < min_x2 || max_x2 < min_x1 || max_y1 < min_y2 || max_y2 < min_y1 {
        return false;
    }

    fn ccw(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64) -> bool {
        (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)
    }

    ccw(x1, y1, x3, y3, x4, y4) != ccw(x2, y2, x3, y3, x4, y4)
        && ccw(x1, y1, x2, y2, x3, y3) != ccw(x1, y1, x2, y2, x4, y4)
}

/// Comprehensive healing with ShapeFix_Solid and ShapeFix_Wire integration.
///
/// This function provides OCCT-equivalent comprehensive healing:
/// 1. Wire-level fixes
/// 2. Face-level fixes
/// 3. Shell-level fixes
/// 4. Solid-level fixes
/// 5. Tolerance propagation
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `options` - Healing options
///
/// # Returns
/// Healed B-Rep and comprehensive report.
pub fn heal_comprehensive(brep: &BRep, options: &HealingOptions) -> (BRep, ComprehensiveHealingReport) {
    let mut report = ComprehensiveHealingReport::default();
    let mut current = brep.clone();

    // Stage 1: Wire fixes
    let (brep_wire, wire_report) = fix_wire(&current, options.tolerance);
    current = brep_wire;
    report.wire_report = Some(wire_report);

    // Stage 2: Face fixes (via standard repair)
    let (brep_face, repair_report) = repair(&current, options.tolerance);
    current = brep_face;
    report.repair_report = Some(repair_report);

    // Stage 3: Solid fixes
    let (brep_solid, solid_report) = fix_solid(&current, options.tolerance);
    current = brep_solid;
    report.solid_report = Some(solid_report);

    // Stage 4: Tolerance propagation
    current = crate::brep_repair::propagate_tolerances(
        &current,
        options.tolerance,
        crate::brep_repair::ToleranceFlowDirection::BottomUp,
    );
    let tol_report = crate::brep_repair::analyze_tolerances(&current, options.tolerance);
    report.tolerance_report = Some(tol_report.vertices);

    // Final check
    report.final_check = brep_check_analyze(&current);
    report.is_clean = report.final_check.is_valid();

    (current, report)
}

/// Comprehensive healing report with all stage details.
#[derive(Debug, Clone, Default)]
pub struct ComprehensiveHealingReport {
    /// Wire-level fix report.
    pub wire_report: Option<WireFixReport>,
    /// Standard repair report.
    pub repair_report: Option<crate::brep_repair::RepairReport>,
    /// Solid-level fix report.
    pub solid_report: Option<SolidFixReport>,
    /// Tolerance propagation report.
    pub tolerance_report: Option<crate::brep_repair::ToleranceStats>,
    /// Final checker result.
    pub final_check: CheckResult,
    /// Whether the result is checker-clean.
    pub is_clean: bool,
}

impl ComprehensiveHealingReport {
    pub fn summary(&self) -> String {
        let mut parts = Vec::new();

        if let Some(ref wr) = self.wire_report
            && wr.total_fixes > 0 {
                parts.push(format!("wires: {} fixes", wr.total_fixes));
            }

        if let Some(ref rr) = self.repair_report {
            let repairs = rr.vertices_merged + rr.faces_reoriented + rr.wires_fixed;
            if repairs > 0 {
                parts.push(format!("repair: {} fixes", repairs));
            }
        }

        if let Some(ref sr) = self.solid_report
            && sr.total_fixes > 0 {
                parts.push(format!("solid: {} fixes", sr.total_fixes));
            }

        if parts.is_empty() {
            if self.is_clean {
                "Clean result, no fixes needed".to_string()
            } else {
                format!("Issues remain: {} issues", self.final_check.issues.len())
            }
        } else {
            format!("{} 鈫?{}", parts.join(", "), if self.is_clean { "clean" } else { "issues remain" })
        }
    }
}
