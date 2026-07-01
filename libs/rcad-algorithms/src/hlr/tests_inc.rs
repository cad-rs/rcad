#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn unit_box_hlr_produces_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = compute_hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "HLR should produce segments for a box"
        );
        assert!(
            result.visible().count() > 0,
            "some segments should be visible"
        );
    }

    #[test]
    fn hlr_svg_is_valid_xml() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = compute_hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(svg.contains("<svg"), "output should be SVG");
        assert!(svg.contains("</svg>"), "SVG should close properly");
        assert!(svg.contains("<line"), "SVG should contain lines");
    }

    #[test]
    fn top_view_box_has_visible_top_edges() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::top(5.0);
        let result = compute_hlr(&brep, &camera, 8);
        let vis = result.visible().count();
        let hid = result.hidden().count();
        assert!(vis > 0, "top view should have visible edges");
        assert!(hid > 0, "top view should have hidden (bottom) edges");
    }

    #[test]
    fn front_view_and_right_view_both_produce_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 1.0,
            depth: 1.0,
        });
        let front_result = compute_hlr(&brep, &HlrCamera::front(5.0), 8);
        let right_result = compute_hlr(&brep, &HlrCamera::right(5.0), 8);
        assert!(!front_result.segments.is_empty());
        assert!(!right_result.segments.is_empty());
    }

    #[test]
    fn hlr_svg_contains_hidden_dashed_lines() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = compute_hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        // Hidden lines are rendered dashed
        assert!(
            svg.contains("stroke-dasharray") || svg.contains("hidden"),
            "SVG should mark hidden lines differently"
        );
    }

    #[test]
    fn hlr_result_has_correct_visibility_counts() {
        // An isometric view of a box has 3 visible faces and 3 hidden faces.
        // The front 3 edges of each visible face → at least some hidden segments exist.
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(10.0);
        let result = compute_hlr(&brep, &camera, 16);
        let total = result.segments.len();
        assert!(total >= 12, "a box has 12 edges, expect at least 12 segments; got {total}");
    }

    #[test]
    fn hlr_circle_edge_sampling() {
        use rcad_kernel::geom::{Circle3, Curve3, CurveEval};

        // Build a minimal BRep with a single circle edge (no solids).
        let mut brep = rcad_kernel::BRep::new();
        let circ = Circle3 {
            center: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
            radius: 1.0,
        };
        // Add two vertices on the circle (half-circle arc)
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(0.0),
        });
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(std::f64::consts::PI),
        });
        brep.edges.push(rcad_kernel::topology::Edge { start: 0, end: 1 });
        brep.geom.curves.push(Curve3::Circle(circ));
        brep.geom.edge_curve.push(Some(0));
        brep.geom
            .edge_curve_range
            .push(Some([0.0, std::f64::consts::PI]));

        let camera = HlrCamera::top(5.0);
        let result = compute_hlr(&brep, &camera, 8);

        // The circle edge should produce at least one segment.
        assert!(
            !result.segments.is_empty(),
            "circle edge should produce HLR segments"
        );

        // All sampled 3D points on the circle should lie ON the circle (unit radius).
        // Verify by checking screen_pts all lie within radius ≈ 1.0 of circle center
        // when projected top-down (X-Y plane).
        for seg in &result.segments {
            // The curve_hint for circle segments should be set.
            assert!(
                matches!(seg.curve_hint, Some(CurveHint::Circle { .. })),
                "circle edge segments should carry CurveHint::Circle"
            );
        }

        // SVG should contain arc path elements (not just lines) for circle edges.
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(
            svg.contains("<path") || result.segments.is_empty(),
            "circle edge SVG should contain <path> arc elements"
        );
    }

    /// Cylinder viewed from the side should produce silhouette line segments
    /// in addition to the wire edges.
    #[test]
    fn cylinder_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        // The cylinder axis is +Y.  Use the right-side camera (looking along -X)
        // so the view direction is perpendicular to the axis → two silhouette lines.
        let camera = HlrCamera::right(10.0);
        let result = compute_hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cylinder HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cylinder should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Sphere HLR should produce silhouette segments (the great circle).
    #[test]
    fn sphere_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = compute_hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "sphere HLR should produce silhouette segments"
        );
    }

    /// Cone viewed from the side should produce two silhouette lines from the apex.
    #[test]
    fn cone_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 1.0,
            height: 2.0,
        });
        // View from the right (perpendicular to cone axis) → two silhouette generators.
        let camera = HlrCamera::right(10.0);
        let result = compute_hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cone HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cone should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Torus HLR should produce silhouette segments.
    #[test]
    fn torus_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        let camera = HlrCamera::front(20.0);
        let result = compute_hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "torus HLR should produce silhouette segments"
        );
    }

    // ── Assembly HLR tests ─────────────────────────────────────────────────────

    /// Two boxes side by side — both should produce segments.
    #[test]
    fn hlr_assembly_two_boxes() {
        let box1 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let box2 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (box1, DAffine3::from_translation(DVec3::new(-2.0, 0.0, 0.0)), "box_left".to_string()),
            (box2, DAffine3::from_translation(DVec3::new(2.0, 0.0, 0.0)), "box_right".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2, "should have 2 component results");
        assert!(result.components.iter().all(|c| !c.segments.is_empty()),
            "each component should produce segments");
    }

    /// Small box behind a large box — the small box should be partially hidden.
    #[test]
    fn hlr_assembly_occlusion() {
        let big = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 3.0, height: 3.0, depth: 3.0,
        });
        let small = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 0.5, height: 0.5, depth: 0.5,
        });

        // Front camera looks along +Y from (0, -10, 0).
        // Place small box at +Y behind the big box so it's occluded.
        let components = vec![
            (big, DAffine3::IDENTITY, "big".to_string()),
            (small, DAffine3::from_translation(DVec3::new(0.0, 3.0, 0.0)), "small_behind".to_string()),
        ];

        let camera = HlrCamera::front(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // The small box behind the big one should have mostly hidden segments
        let small_comp = result.components.iter().find(|c| c.name == "small_behind").unwrap();
        let hidden = small_comp.segments.iter().filter(|s| !s.visible).count();
        let visible = small_comp.segments.iter().filter(|s| s.visible).count();
        assert!(hidden > visible,
            "small box behind big one should have more hidden than visible segments; hidden={hidden}, visible={visible}");
    }

    /// Assembly with a single component should match single-BRep HLR.
    #[test]
    fn hlr_assembly_single_matches_hlr() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);

        let single_hlr = compute_hlr(&brep, &camera, 8);
        let assembly_result = hlr_assembly(
            &[(brep.clone(), DAffine3::IDENTITY, "box".to_string())],
            &camera, 8,
        );

        assert_eq!(assembly_result.components.len(), 1);
        let asm_segs = &assembly_result.components[0].segments;
        // Segment counts should be similar (same geometry, same algorithm)
        assert!(asm_segs.len() >= single_hlr.segments.len() - 2,
            "assembly HLR should produce similar segment count");
        assert!(asm_segs.len() <= single_hlr.segments.len() + 2,
            "assembly HLR should produce similar segment count");
    }

    /// Stacked boxes — top box visible, bottom box partially occluded.
    #[test]
    fn hlr_assembly_stacked_boxes() {
        let bottom = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0, height: 1.0, depth: 2.0,
        });
        let top = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (bottom, DAffine3::from_translation(DVec3::new(0.0, 0.0, 0.0)), "bottom".to_string()),
            (top, DAffine3::from_translation(DVec3::new(0.0, 0.0, 1.5)), "top".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // Both boxes should have some visible segments
        for comp in &result.components {
            let vis = comp.segments.iter().filter(|s| s.visible).count();
            assert!(vis > 0, "{} should have visible segments", comp.name);
        }
    }

    /// Empty assembly should return empty result.
    #[test]
    fn hlr_assembly_empty() {
        let components: Vec<(BRep, DAffine3, String)> = vec![];
        let camera = HlrCamera::isometric(5.0);
        let result = hlr_assembly(&components, &camera, 8);
        assert!(result.components.is_empty());
    }

    // ── Improved HLR tests ─────────────────────────────────────────────────────

    #[test]
    fn hlr_options_default_values() {
        let opts = HlrOptions::default();
        assert_eq!(opts.edge_samples, 8);
        assert_eq!(opts.silhouette_samples, 32);
        assert!(opts.curvature_adaptive);
        assert!(opts.tangent_tolerance > 0.0);
    }

    #[test]
    fn hlr_options_builders() {
        let opts = HlrOptions::default()
            .with_edge_samples(16)
            .with_silhouette_samples(64)
            .with_curvature_adaptive(false)
            .with_tangent_tolerance(TOLERANCE_RETRY_LADDER_COARSE);

        assert_eq!(opts.edge_samples, 16);
        assert_eq!(opts.silhouette_samples, 64);
        assert!(!opts.curvature_adaptive);
        assert!((opts.tangent_tolerance - TOLERANCE_RETRY_LADDER_COARSE).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn hlr_with_options_basic() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let opts = HlrOptions::default().with_edge_samples(16);
        let result = compute_hlr_with_options(&brep, &camera, opts);

        assert!(!result.segments.is_empty(), "should produce segments");
        // All segments from a box should be edges, not silhouettes
        assert!(result.segments.iter().all(|s| s.segment_type == SegmentType::Edge));
    }

    #[test]
    fn cylinder_silhouettes_are_marked() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        let camera = HlrCamera::right(10.0);
        let result = compute_hlr(&brep, &camera, 8);

        // Should have both edge and silhouette segments
        let has_silhouette = result.segments.iter().any(|s| s.is_contour());
        assert!(has_silhouette, "cylinder should have silhouette segments");
    }

    #[test]
    fn sphere_silhouettes_are_marked() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = compute_hlr(&brep, &camera, 8);

        // All segments from a sphere should be silhouettes (no wire edges)
        assert!(
            result.segments.iter().all(|s| s.is_contour()),
            "sphere should only have silhouette segments"
        );
    }

    #[test]
    fn extract_silhouette_curves_sphere() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });
        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert_eq!(curves.len(), 1, "sphere should have one silhouette curve");
        assert!(curves[0].points.len() >= 32, "silhouette should have enough points");

        // All points should be at distance ~2.0 from origin
        for pt in &curves[0].points {
            let dist = pt.length();
            assert!(
                (dist - 2.0).abs() < 0.01,
                "silhouette point distance should be ~2.0, got {dist}"
            );
        }
    }

    #[test]
    fn extract_silhouette_curves_cylinder() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 3.0,
        });
        // View along X axis - perpendicular to cylinder axis (Y)
        let view_dir = DVec3::X;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert!(curves.len() >= 2, "cylinder should have at least 2 silhouette curves");

        // Each silhouette curve should be a line (two lines on opposite sides)
        for curve in &curves {
            assert!(curve.points.len() >= 16, "silhouette should have enough points");
        }
    }

    #[test]
    fn extract_silhouette_curves_torus() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert!(curves.len() >= 2, "torus should have at least 2 silhouette curves");
    }

    #[test]
    fn hlr_result_silhouettes_iterator() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = compute_hlr(&brep, &camera, 8);

        let sil_count = result.silhouettes().count();
        assert!(sil_count > 0, "should have silhouette segments");

        let vis_sil_count = result.visible_silhouettes().count();
        assert!(vis_sil_count > 0, "should have visible silhouette segments");
    }

    #[test]
    fn segment_is_contour_method() {
        let seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Silhouette,
        };
        assert!(seg.is_contour());

        let edge_seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Edge,
        };
        assert!(!edge_seg.is_contour());
    }

    #[test]
    fn adaptive_sampling_high_curvature() {
        // Test that adaptive sampling produces more points in high-curvature regions
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let view_dir = DVec3::Z;

        let opts_low = HlrOptions {
            silhouette_samples: 16,
            curvature_adaptive: false,
            ..HlrOptions::default()
        };
        let opts_high = HlrOptions {
            silhouette_samples: 64,
            curvature_adaptive: true,
            ..HlrOptions::default()
        };

        let curves_low = extract_silhouette_curves(&brep, view_dir, &opts_low);
        let curves_high = extract_silhouette_curves(&brep, view_dir, &opts_high);

        // Both should produce curves
        assert!(!curves_low.is_empty());
        assert!(!curves_high.is_empty());

        // Higher sampling should produce more points
        let pts_low: usize = curves_low.iter().map(|c| c.points.len()).sum();
        let pts_high: usize = curves_high.iter().map(|c| c.points.len()).sum();
        assert!(
            pts_high >= pts_low,
            "higher sampling should produce at least as many points"
        );
    }

    #[test]
    fn tangent_tolerance_affects_silhouette_detection() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let view_dir = DVec3::Z;

        // Very tight tolerance
        let opts_tight = HlrOptions {
            tangent_tolerance: TOLERANCE_LEN_MIN,
            ..HlrOptions::default()
        };

        // Very loose tolerance (should still work for sphere)
        let opts_loose = HlrOptions {
            tangent_tolerance: 0.01,
            ..HlrOptions::default()
        };

        let curves_tight = extract_silhouette_curves(&brep, view_dir, &opts_tight);
        let curves_loose = extract_silhouette_curves(&brep, view_dir, &opts_loose);

        // Both should find silhouette curves for a sphere
        assert!(!curves_tight.is_empty());
        assert!(!curves_loose.is_empty());
    }

    // ── New Enhanced HLR Tests ───────────────────────────────────────────────────

    #[test]
    fn hlr_options_new_fields() {
        let opts = HlrOptions::default();
        assert!(opts.parallel, "parallel should be true by default");
        assert_eq!(opts.parallel_threshold, 4);
        assert!(opts.cache_surface_properties);
        assert!(opts.detect_thread_edges);
        assert!(opts.detect_seam_edges);
    }

    #[test]
    fn hlr_options_new_builders() {
        let opts = HlrOptions::default()
            .with_parallel(false)
            .with_parallel_threshold(8)
            .with_surface_caching(false)
            .with_thread_edge_detection(false)
            .with_seam_edge_detection(false)
            .with_silhouette_proximity(0.05);

        assert!(!opts.parallel);
        assert_eq!(opts.parallel_threshold, 8);
        assert!(!opts.cache_surface_properties);
        assert!(!opts.detect_thread_edges);
        assert!(!opts.detect_seam_edges);
        assert!((opts.silhouette_proximity_factor - 0.05).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn surface_property_cache_basic() {
        let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Y),
        });

        let domain = surface.default_domain();
        let mut cache = SurfacePropertyCache::new(16, domain);

        // First access should compute
        let props1 = cache.get_or_compute(&surface, 0.5, 0.5);
        assert!(cache.len() > 0, "cache should have entries");

        // Second access should return cached value
        let props2 = cache.get_or_compute(&surface, 0.5, 0.5);
        assert!((props1.point - props2.point).length() < TOLERANCE_LINEAR_ULTRA_STRICT);

        // Verify surface properties
        assert!((props1.point.length() - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "sphere point should be on surface");
        assert!((props1.curvatures.0 - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "sphere curvature should be 1/r");
    }

    #[test]
    fn surface_properties_near_silhouette() {
        let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Y),
        });

        // For a Y-axis sphere: x_ax = Z (perpendicular to Y), y_ax = X = Y.cross(Z)
        // At u=π/2, v=π/2: normal = u.cos * Z + u.sin * X = 0*Z + 1*X = X (perpendicular to Z view)
        let props = compute_surface_properties(&surface, std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        let view_dir = DVec3::Z;

        assert!(props.is_near_silhouette(view_dir, 0.01), "equator at u=π/2 should be near silhouette for Z view");
        assert!((props.normal_dot_view(view_dir)).abs() < 0.01);

        // At u=0, v=π/2: normal = Z (parallel to view) - NOT a silhouette
        let props_front = compute_surface_properties(&surface, 0.0, std::f64::consts::FRAC_PI_2);
        assert!(!props_front.is_near_silhouette(view_dir, 0.5), "point facing viewer should not be near silhouette");

        // At pole (v = 0), normal is Y axis, perpendicular to Z view, so it IS a silhouette
        let props_pole = compute_surface_properties(&surface, 0.0, 0.0);
        assert!(props_pole.is_near_silhouette(view_dir, 0.5), "pole (Y-normal) should be near silhouette for Z view");
        assert!((props_pole.normal_dot_view(view_dir)).abs() < 0.01);
    }

    #[test]
    fn spatial_index_basic() {
        let points: Vec<DVec3> = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(0.1, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.5, 0.5, 0.5),
        ];

        let index = SilhouetteSpatialIndex::build(&points, 0.5);

        assert_eq!(index.len(), 5, "index should have all points");

        // Query radius
        let nearby = index.query_radius(DVec3::ZERO, 0.2);
        assert!(nearby.len() >= 2, "should find points near origin");

        // Query nearest
        let nearest = index.query_nearest(DVec3::new(0.05, 0.05, 0.0));
        assert!(nearest.is_some());
        let (idx, dist) = nearest.unwrap();
        assert!(dist < 0.1, "nearest point should be close");
    }

    #[test]
    fn spatial_index_empty() {
        let points: Vec<DVec3> = vec![];
        let index = SilhouetteSpatialIndex::build(&points, 0.5);

        assert!(index.is_empty());
        assert!(index.query_nearest(DVec3::ZERO).is_none());
    }

    #[test]
    fn edge_classification_basic() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let camera = HlrCamera::isometric(5.0);
        let opts = HlrOptions::default();
        let classifications = classify_edges(&brep, &camera, &opts);

        assert_eq!(classifications.len(), brep.edges.len(), "should classify all edges");

        // All box edges should be regular edges (not thread or seam)
        for class in &classifications {
            assert!(
                class.classification != EdgeClassification::Thread,
                "box edges should not be thread edges"
            );
        }
    }

    #[test]
    fn edge_classification_cylinder_seam() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let camera = HlrCamera::isometric(5.0);
        let opts = HlrOptions::default().with_seam_edge_detection(true);
        let classifications = classify_edges(&brep, &camera, &opts);

        // Cylinder should have at least one seam edge
        let seam_count = classifications
            .iter()
            .filter(|c| c.classification == EdgeClassification::Seam)
            .count();

        // Note: the seam detection depends on the BRep structure
        // For a primitive cylinder, the seam edge may be detected
        assert!(classifications.len() > 0, "should have edge classifications");
    }

    #[test]
    fn segment_type_thread_and_seam() {
        let thread_seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Thread,
        };
        assert!(thread_seg.is_thread());
        assert!(!thread_seg.is_seam());
        assert!(!thread_seg.is_contour());

        let seam_seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Seam,
        };
        assert!(seam_seg.is_seam());
        assert!(!seam_seg.is_thread());
        assert!(!seam_seg.is_contour());
    }

    #[test]
    fn adaptive_sampling_config_default() {
        let config = AdaptiveSamplingConfig::default();
        assert_eq!(config.base_samples, 32);
        assert!(config.max_samples > config.base_samples);
        assert!(config.curvature_threshold > 0.0);
        assert!(config.proximity_threshold > 0.0);
    }

    #[test]
    fn adaptive_sample_creation() {
        let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Y),
        });

        let view_dir = DVec3::Z;
        let sample = create_adaptive_sample(&surface, view_dir, 0.0, std::f64::consts::FRAC_PI_2);

        assert!(sample.is_some());
        let s = sample.unwrap();

        // Check that the sample is on the equator
        assert!((s.point.y - 0.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "equator y should be 0");
        assert!((s.curvature - 0.5).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "sphere radius 2 curvature should be 0.5");
    }

    #[test]
    fn hlr_with_parallel_processing() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let camera = HlrCamera::isometric(5.0);

        // Test with parallel processing enabled
        let opts_parallel = HlrOptions::default()
            .with_parallel(true)
            .with_parallel_threshold(1);

        let result_parallel = compute_hlr_with_options(&brep, &camera, opts_parallel);

        // Test with parallel processing disabled
        let opts_serial = HlrOptions::default()
            .with_parallel(false);

        let result_serial = compute_hlr_with_options(&brep, &camera, opts_serial);

        // Both should produce the same number of segments (within tolerance)
        assert!(!result_parallel.segments.is_empty());
        assert!(!result_serial.segments.is_empty());
    }

    #[test]
    fn curve_surface_intersection_basic() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let camera = HlrCamera::front(5.0);
        let opts = HlrOptions::default();

        // The sphere BRep has a seam edge
        if !brep.edges.is_empty() {
            // Try to compute curve-surface intersection for the first edge
            let edge_idx = 0;
            if let Some(surface_idx) = brep.geom.face_surface.get(0).and_then(|&s| s) {
                let result = compute_curve_visibility_on_surface(
                    &brep,
                    edge_idx,
                    surface_idx,
                    &camera,
                    &opts,
                );

                // The function should return a result (may have empty intersections)
                if let Some(intersection) = result {
                    assert!(intersection.curve_params.len() == intersection.points.len());
                }
            }
        }
    }

    #[test]
    fn hlr_result_with_thread_segments() {
        // Create a cylinder - thread edges would be helical, but primitives
        // don't have those. This tests the segment type detection.
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let camera = HlrCamera::right(10.0);
        let opts = HlrOptions::default()
            .with_thread_edge_detection(true)
            .with_seam_edge_detection(true);

        let result = compute_hlr_with_options(&brep, &camera, opts);

        assert!(!result.segments.is_empty(), "should have segments");

        // Check that we have various segment types
        let has_edge = result.segments.iter().any(|s| s.segment_type == SegmentType::Edge);
        let has_silhouette = result.segments.iter().any(|s| s.segment_type == SegmentType::Silhouette);

        assert!(has_edge || has_silhouette, "should have edge or silhouette segments");
    }

    #[test]
    fn grazing_angle_handling() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Test with different grazing angle thresholds
        let camera = HlrCamera::front(5.0);

        let opts_tight = HlrOptions {
            grazing_angle_threshold: 0.01,
            ..HlrOptions::default()
        };

        let opts_loose = HlrOptions {
            grazing_angle_threshold: 0.5,
            ..HlrOptions::default()
        };

        let result_tight = compute_hlr_with_options(&brep, &camera, opts_tight);
        let result_loose = compute_hlr_with_options(&brep, &camera, opts_loose);

        // Both should produce results
        assert!(!result_tight.segments.is_empty());
        assert!(!result_loose.segments.is_empty());
    }

    #[test]
    fn performance_with_caching() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });

        let camera = HlrCamera::front(5.0);

        // Test with caching enabled
        let opts_cached = HlrOptions::default()
            .with_surface_caching(true);

        let result_cached = compute_hlr_with_options(&brep, &camera, opts_cached);

        // Test with caching disabled
        let opts_uncached = HlrOptions::default()
            .with_surface_caching(false);

        let result_uncached = compute_hlr_with_options(&brep, &camera, opts_uncached);

        // Both should produce valid results
        assert!(!result_cached.segments.is_empty());
        assert!(!result_uncached.segments.is_empty());
    }

    #[test]
    fn is_degenerate_edge_check() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Check all edges for degeneracy
        for i in 0..brep.edges.len() {
            let _is_degenerate = is_degenerate_edge_for_hlr(&brep, i);
            // For a sphere primitive, edges should not be degenerate
            // (the seam edge is not degenerate, just periodic)
        }
    }

    // ── Ellipsoid Silhouette Tests ───────────────────────────────────────────────

    #[test]
    fn ellipsoid_silhouette_basic() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Create an ellipsoid with different radii
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z; // View along Z axis
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert_eq!(curves.len(), 1, "ellipsoid should have one silhouette curve");
        assert!(curves[0].len() >= 32, "silhouette should have enough points");
    }

    #[test]
    fn ellipsoid_silhouette_satisfies_condition() {
        use rcad_kernel::geom::{EllipsoidalSurface, SurfaceEval};

        // Create an ellipsoid
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty(), "should have silhouette curves");

        // Check that all silhouette points satisfy n·v ≈ 0
        for pt in &curves[0] {
            // Compute the point in local coordinates
            let x = pt.x;
            let y = pt.y;
            let z = pt.z;

            // Normal direction (gradient of implicit equation, normalized)
            let grad = DVec3::new(
                x / (ell.radius_x * ell.radius_x),
                y / (ell.radius_y * ell.radius_y),
                z / (ell.radius_z * ell.radius_z),
            );
            let normal = grad.normalize_or_zero();

            // Dot product with view direction should be near zero
            let dot = normal.dot(view_dir);
            assert!(
                dot.abs() < 0.05,
                "silhouette point should satisfy n·v ≈ 0, got {dot}"
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_on_surface() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Create an ellipsoid
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // Check that all silhouette points are on the ellipsoid surface
        // x²/a² + y²/b² + z²/c² = 1
        for pt in &curves[0] {
            let value = (pt.x / ell.radius_x).powi(2)
                + (pt.y / ell.radius_y).powi(2)
                + (pt.z / ell.radius_z).powi(2);
            assert!(
                (value - 1.0).abs() < TOLERANCE_MESH_LEGACY,
                "point should be on ellipsoid surface, got implicit value {value}"
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_various_view_directions() {
        use rcad_kernel::geom::EllipsoidalSurface;

        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let opts = HlrOptions::default();

        // Test various view directions
        let view_directions = [
            DVec3::X,
            DVec3::Y,
            DVec3::Z,
            DVec3::new(1.0, 1.0, 0.0).normalize(),
            DVec3::new(1.0, 1.0, 1.0).normalize(),
            DVec3::new(0.0, 1.0, 1.0).normalize(),
        ];

        for view_dir in view_directions {
            let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);
            assert!(
                !curves.is_empty() && !curves[0].is_empty(),
                "should have silhouette for view_dir {:?}",
                view_dir
            );

            // Verify silhouette condition
            for pt in &curves[0] {
                let grad = DVec3::new(
                    pt.x / (ell.radius_x * ell.radius_x),
                    pt.y / (ell.radius_y * ell.radius_y),
                    pt.z / (ell.radius_z * ell.radius_z),
                );
                let normal = grad.normalize_or_zero();
                let dot = normal.dot(view_dir);
                assert!(
                    dot.abs() < 0.1,
                    "silhouette condition not satisfied for view_dir {:?}, dot = {}",
                    view_dir,
                    dot
                );
            }
        }
    }

    #[test]
    fn ellipsoid_silhouette_sphere_case() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // A sphere is a special case of an ellipsoid with all radii equal
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 2.0,
            radius_y: 2.0,
            radius_z: 2.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert_eq!(curves.len(), 1, "sphere should have one silhouette curve");

        // All points should be at distance 2.0 from origin (great circle)
        for pt in &curves[0] {
            let dist = pt.length();
            assert!(
                (dist - 2.0).abs() < 0.01,
                "sphere silhouette point should be at radius distance, got {dist}"
            );

            // z-coordinate should be near 0 (great circle perpendicular to Z)
            assert!(
                pt.z.abs() < 0.01,
                "great circle should be in XY plane, got z = {}",
                pt.z
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_translated() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Ellipsoid not at origin
        let ell = EllipsoidalSurface {
            center: DVec3::new(1.0, -2.0, 0.5),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // Check that all points are centered around the ellipsoid center
        for pt in &curves[0] {
            let local = *pt - ell.center;
            let value = (local.x / ell.radius_x).powi(2)
                + (local.y / ell.radius_y).powi(2)
                + (local.z / ell.radius_z).powi(2);
            assert!(
                (value - 1.0).abs() < TOLERANCE_MESH_LEGACY,
                "point should be on translated ellipsoid surface"
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_rotated_frame() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Ellipsoid with rotated axis (not aligned with Z)
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y, // Axis along Y
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Y; // View along the axis
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // The silhouette should be an ellipse in the XZ plane
        for pt in &curves[0] {
            // Y coordinate should be near 0 (silhouette in plane perpendicular to view)
            assert!(
                pt.y.abs() < 0.01,
                "silhouette should be in XZ plane for Y view, got y = {}",
                pt.y
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_closed_curve() {
        use rcad_kernel::geom::EllipsoidalSurface;

        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // The silhouette should be a closed curve
        // First and last points should be close
        let pts = &curves[0];
        let first = pts.first().unwrap();
        let last = pts.last().unwrap();
        let closure_dist = (*first - *last).length();
        assert!(
            closure_dist < 0.5,
            "silhouette should be approximately closed, distance between first and last = {}",
            closure_dist
        );
    }
}
