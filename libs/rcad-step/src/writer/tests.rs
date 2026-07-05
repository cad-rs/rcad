#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        StepDatum, StepDimensionalLocation, StepDimensionalSize, StepGeometricTolerance,
        StepGeometricToleranceWithDatumReference,
        StepPropertyDefinitionRepr, StepReader, ToleranceZonePosition, ToleranceZoneShape,
    };
    use glam::DVec3;
    use rcad_modeling::make_box_brep;
    use std::io::Cursor;
    const HFSS_STEP: &str = include_str!("../../../assets/hfss.step");

    /// Convert topods::BRep to FlatBRep for field access in tests.
    fn to_flat(t: &rcad_kernel::topods::BRep) -> crate::brep_flat::FlatBRep {
        crate::brep_flat::FlatBRep::from_topods(t)
    }

    #[test]
    fn exports_full_box_and_reimports() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0)
            .expect("test box should be valid");
        // Default export path favors standard solid-only representation.
        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("exported STEP should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);
        assert!(!reparsed.edges.is_empty());
        assert!(!reparsed.solids.is_empty());
        assert!(step.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"));
        assert!(step.contains("MANIFOLD_SOLID_BREP"));
        assert!(step.contains("CLOSED_SHELL"));
        // Full solid export should not include wireframe overlays by default.
        assert!(!step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"));
        assert!(!step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("SHAPE_DEFINITION_REPRESENTATION"));
    }

    #[test]
    fn exports_general_properties_when_provided() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let props = vec![
            StepGeneralProperty {
                name: "PartNumber".to_string(),
                description: Some("PN-001".to_string()),
            },
            StepGeneralProperty {
                name: "Revision".to_string(),
                description: Some("A".to_string()),
            },
        ];
        let step = StepWriter::write_string_with_properties(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &props,
            StepProtocol::Ap242,
        );

        assert!(step.contains("GENERAL_PROPERTY('PartNumber','PN-001',$)"));
        assert!(step.contains("GENERAL_PROPERTY('Revision','A',$)"));
        assert!(step.contains("PROPERTY_DEFINITION('PartNumber','PN-001',#"));
        assert!(step.contains("PROPERTY_DEFINITION('Revision','A',#"));
    }

    #[test]
    fn exports_ap242_metadata_entities_when_provided() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            property_definition_representations: vec![StepPropertyDefinitionRepr {
                entity_id: 0,
                property_definition_id: Some(10),
                representation_id: Some(20),
            }],
            dimensional_locations: vec![StepDimensionalLocation {
                entity_id: 0,
                name: Some("d_loc".into()),
                description: Some("desc".into()),
                from_entity_id: Some(30),
                to_entity_id: Some(31),
            }],
            dimensional_sizes: vec![StepDimensionalSize {
                entity_id: 0,
                name: Some("d_size".into()),
                description: None,
                shape_aspect_id: Some(40),
            }],
            geometric_tolerances: vec![StepGeometricTolerance {
                entity_id: 0,
                name: Some("flatness".into()),
                description: Some("gtol".into()),
                value_entity_id: Some(50),
                shape_aspect_id: Some(60),
            }],
            geometric_tolerances_with_datum_references: vec![
                StepGeometricToleranceWithDatumReference {
                    entity_id: 0,
                    name: Some("position".into()),
                    description: Some("gtol_datum".into()),
                    value_entity_id: Some(51),
                    shape_aspect_id: Some(61),
                    datum_system_id: Some(71),
                },
            ],
            datums: vec![StepDatum {
                entity_id: 0,
                name: Some("A".into()),
                description: Some("primary".into()),
                shape_aspect_id: Some(70),
            }],
            datum_systems: vec![StepDatumSystem {
                entity_id: 0,
                name: Some("A_SYS".into()),
                description: Some("primary_system".into()),
                datum_ids: vec![70],
            }],
            kinematic_pairs: vec![StepKinematicPair {
                entity_id: 0,
                entity_type: "REVOLUTE_PAIR".into(),
                name: Some("hinge".into()),
                description: Some("joint".into()),
                related_entity_ids: vec![81, 82, 83],
            }],
            // GDT extended fields
            dimensional_tolerances: vec![],
            tolerance_values: vec![],
            position_tolerances: vec![],
            orientation_tolerances: vec![],
            form_tolerances: vec![],
            runout_tolerances: vec![],
            profile_tolerances: vec![],
            datum_reference_frames: vec![],
            datum_targets: vec![],
            tolerance_zone_definitions_enhanced: vec![],
            // View and annotation fields
            views: vec![],
            cameras: vec![],
            view_volumes: vec![],
            notes: vec![],
            annotation_planes: vec![],
            annotation_occurrences: vec![],
            dimension_curves: vec![],
            terminator_symbols: vec![],
            datum_feature_callouts: vec![],
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );

        assert!(step.contains("PROPERTY_DEFINITION_REPRESENTATION(#10,#20)"));
        assert!(step.contains("DIMENSIONAL_LOCATION('d_loc','desc',#30,#31)"));
        assert!(step.contains("DIMENSIONAL_SIZE('d_size',$,#40)"));
        assert!(step.contains("GEOMETRIC_TOLERANCE('flatness','gtol',#50,#60)"));
        assert!(step.contains(
            "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE('position','gtol_datum',#51,#61,#71)"
        ));
        assert!(step.contains("DATUM('A','primary',#70)"));
        assert!(step.contains("DATUM_SYSTEM('A_SYS','primary_system',(#70))"));
        assert!(step.contains("REVOLUTE_PAIR('hinge','joint',#81,#82,#83)"));
    }

    #[test]
    fn ap242_metadata_write_read_roundtrip_smoke() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            property_definition_representations: vec![StepPropertyDefinitionRepr {
                entity_id: 0,
                property_definition_id: Some(11),
                representation_id: Some(22),
            }],
            dimensional_locations: vec![StepDimensionalLocation {
                entity_id: 0,
                name: Some("L1".into()),
                description: Some("loc".into()),
                from_entity_id: Some(33),
                to_entity_id: Some(44),
            }],
            dimensional_sizes: vec![StepDimensionalSize {
                entity_id: 0,
                name: Some("S1".into()),
                description: Some("size".into()),
                shape_aspect_id: Some(55),
            }],
            geometric_tolerances: vec![StepGeometricTolerance {
                entity_id: 0,
                name: Some("parallelism".into()),
                description: Some("tol".into()),
                value_entity_id: Some(66),
                shape_aspect_id: Some(77),
            }],
            geometric_tolerances_with_datum_references: vec![
                StepGeometricToleranceWithDatumReference {
                    entity_id: 0,
                    name: Some("perpendicularity".into()),
                    description: Some("tol_datum".into()),
                    value_entity_id: Some(67),
                    shape_aspect_id: Some(78),
                    datum_system_id: Some(89),
                },
            ],
            datums: vec![StepDatum {
                entity_id: 0,
                name: Some("B".into()),
                description: Some("secondary".into()),
                shape_aspect_id: Some(88),
            }],
            datum_systems: vec![StepDatumSystem {
                entity_id: 0,
                name: Some("B_SYS".into()),
                description: Some("secondary_system".into()),
                datum_ids: vec![88],
            }],
            kinematic_pairs: vec![StepKinematicPair {
                entity_id: 0,
                entity_type: "PRISMATIC_PAIR".into(),
                name: Some("slider".into()),
                description: Some("guide".into()),
                related_entity_ids: vec![90, 91],
            }],
            // GDT extended fields
            dimensional_tolerances: vec![StepDimensionalTolerance {
                entity_id: 0,
                name: Some("diam_tol".into()),
                description: Some("diameter tolerance".into()),
                dimensional_characteristic_id: Some(100),
                upper_tolerance: Some(0.05),
                lower_tolerance: Some(-0.05),
                unit: Some("mm".into()),
            }],
            tolerance_values: vec![StepToleranceValue {
                entity_id: 0,
                name: Some("tol_val".into()),
                value: 0.025,
                unit: Some("mm".into()),
            }],
            position_tolerances: vec![StepPositionTolerance {
                entity_id: 0,
                name: Some("pos_tol".into()),
                description: Some("positional tolerance".into()),
                value_entity_id: Some(101),
                shape_aspect_id: Some(102),
                datum_system_id: Some(103),
                projected: false,
                projected_height: None,
            }],
            orientation_tolerances: vec![StepOrientationTolerance {
                entity_id: 0,
                name: Some("ang_tol".into()),
                description: Some("angularity".into()),
                value_entity_id: Some(104),
                shape_aspect_id: Some(105),
                datum_system_id: Some(106),
                orientation_type: OrientationToleranceType::Angularity,
            }],
            form_tolerances: vec![StepFormTolerance {
                entity_id: 0,
                name: Some("flat_tol".into()),
                description: Some("flatness".into()),
                value_entity_id: Some(107),
                shape_aspect_id: Some(108),
                form_type: FormToleranceType::Flatness,
            }],
            runout_tolerances: vec![StepRunoutTolerance {
                entity_id: 0,
                name: Some("cr_tol".into()),
                description: Some("circular runout".into()),
                value_entity_id: Some(109),
                shape_aspect_id: Some(110),
                datum_system_id: Some(111),
                runout_type: RunoutToleranceType::CircularRunout,
            }],
            profile_tolerances: vec![StepProfileTolerance {
                entity_id: 0,
                name: Some("lin_tol".into()),
                description: Some("profile of a line".into()),
                value_entity_id: Some(112),
                shape_aspect_id: Some(113),
                datum_system_id: None,
                profile_type: ProfileToleranceType::ProfileOfALine,
            }],
            datum_reference_frames: vec![StepDatumReferenceFrame {
                entity_id: 0,
                name: Some("DRF1".into()),
                description: Some("datum reference frame".into()),
                datum_system_ids: vec![88],
            }],
            datum_targets: vec![StepDatumTarget {
                entity_id: 0,
                name: Some("A1".into()),
                description: Some("datum target".into()),
                target_identifier: Some("A1".into()),
                datum_id: Some(88),
                target_type: DatumTargetType::Point,
                shape_aspect_id: Some(114),
            }],
            tolerance_zone_definitions_enhanced: vec![StepToleranceZoneDefinitionEnhanced {
                entity_id: 0,
                name: Some("cylindrical".into()),
                description: Some("symmetric".into()),
                tolerance_zone_id: Some(115),
                shape_aspect_id: Some(116),
                zone_shape: ToleranceZoneShape::Cylindrical,
                zone_position: ToleranceZonePosition::Symmetric,
                defining_shape_aspect_id: None,
            }],
            // View and annotation fields
            views: vec![],
            cameras: vec![],
            view_volumes: vec![],
            notes: vec![],
            annotation_planes: vec![],
            annotation_occurrences: vec![],
            dimension_curves: vec![],
            terminator_symbols: vec![],
            datum_feature_callouts: vec![],
        };

        let step = StepWriter::write_string_with_ap242_metadata(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );

        let (_parsed_brep, doc_meta) =
            StepReader::parse_string_with_metadata(&step).expect("AP242 metadata STEP should parse");
        assert_eq!(doc_meta.property_definition_representations.len(), 1);
        assert_eq!(doc_meta.dimensional_locations.len(), 1);
        assert_eq!(doc_meta.dimensional_sizes.len(), 1);
        assert_eq!(doc_meta.geometric_tolerances.len(), 1);
        assert_eq!(doc_meta.geometric_tolerances_with_datum_references.len(), 1);
        assert_eq!(doc_meta.datums.len(), 1);
        assert_eq!(doc_meta.datum_systems.len(), 1);
        assert_eq!(doc_meta.kinematic_pairs.len(), 1);
        // GDT extended assertions
        assert_eq!(doc_meta.dimensional_tolerances.len(), 1);
        assert_eq!(doc_meta.tolerance_values.len(), 1);
        assert_eq!(doc_meta.position_tolerances.len(), 1);
        assert_eq!(doc_meta.orientation_tolerances.len(), 1);
        assert_eq!(doc_meta.form_tolerances.len(), 1);
        assert_eq!(doc_meta.runout_tolerances.len(), 1);
        assert_eq!(doc_meta.profile_tolerances.len(), 1);
        assert_eq!(doc_meta.datum_reference_frames.len(), 1);
        assert_eq!(doc_meta.datum_targets.len(), 1);
        assert_eq!(doc_meta.tolerance_zone_definitions_enhanced.len(), 1);
    }

    #[test]
    fn exports_selected_edges_without_faces() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[0, 1],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("edge-only export should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);
        assert!(reparsed.solids.is_empty());
        assert_eq!(reparsed.edges.len(), 2);
    }

    #[test]
    fn standalone_reversed_range_exports_false_sense_trimmed_curve() {
        let mut brep = BRep::new();
        brep.vertices = vec![
            Vertex {
                point: DVec3::new(1.0, 0.0, 0.0),
            },
            Vertex {
                point: DVec3::new(0.0, 1.0, 0.0),
            },
        ];
        brep.edges = vec![Edge { start: 0, end: 1 }];
        brep.geom.curves.push(Curve3::Circle(rcad_kernel::geom::Circle3::new(DVec3::ZERO, DVec3::Z, 1.0,
        )));
        brep.geom.edge_curve = vec![Some(0)];
        brep.geom.edge_curve_range = vec![Some([135.0, 0.0])];

        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[0],
            },
        );
        assert!(step.contains("TRIMMED_CURVE("));
        assert!(step.contains(",.F.,.PARAMETER.)"));
    }

    #[test]
    fn standalone_circle_uses_degree_range_hint_for_major_arc_sweep() {
        let mut brep = BRep::new();
        brep.vertices = vec![
            Vertex {
                point: DVec3::new(1.0, 0.0, 0.0),
            },
            Vertex {
                point: DVec3::new(0.0, 1.0, 0.0),
            },
        ];
        brep.edges = vec![Edge { start: 0, end: 1 }];
        brep.geom.curves.push(Curve3::Circle(rcad_kernel::geom::Circle3::new(DVec3::ZERO, DVec3::Z, 1.0,
        )));
        brep.geom.edge_curve = vec![Some(0)];
        // 270-degree sweep hint should choose the major arc.
        brep.geom.edge_curve_range = vec![Some([0.0, 270.0])];

        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[0],
            },
        );
        assert!(step.contains("TRIMMED_CURVE("));
        assert!(step.contains("PARAMETER_VALUE(270.000"));
    }

    #[test]
    fn stream_write_then_stream_read_roundtrip() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0)
            .expect("test box should be valid");

        let mut buf = Vec::<u8>::new();
        StepWriter::write_to(
            &mut buf,
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        )
        .expect("stream write should succeed");

        let reparsed =
            StepReader::parse_reader(Cursor::new(buf)).expect("stream read should parse");
        let reparsed = to_flat(&reparsed);
        assert!(!reparsed.edges.is_empty());
        assert!(!reparsed.solids.is_empty());
    }

    #[test]
    fn exports_selected_faces_via_shell_based_surface_model() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[0],
                selected_edges: &[],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("selected-face export should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);
        assert!(!reparsed.solids.is_empty());
        assert!(step.contains("OPEN_SHELL"));
        assert!(step.contains("SHELL_BASED_SURFACE_MODEL"));
        assert!(step.contains("MANIFOLD_SURFACE_SHAPE_REPRESENTATION"));
    }

    #[test]
    fn selected_faces_also_export_boundary_wire_entities() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string_with_options(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[0],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                ..Default::default()
            },
        );

        assert!(
            step.contains("GEOMETRIC_CURVE_SET"),
            "selected-face export should include boundary 1D entities"
        );
        assert!(
            step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"),
            "selected-face export should include wireframe representation"
        );
    }

    #[test]
    fn exports_analytic_surfaces_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse"); let brep_t = brep; let brep = to_flat(&brep_t);
        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        // All analytic surfaces should now be exported properly as ADVANCED_FACE
        // with their respective surface types, including seam faces on
        // spheres and cones.
        assert!(step.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"));
        assert!(step.contains("MANIFOLD_SURFACE_SHAPE_REPRESENTATION"));
        assert!(step.contains("SPHERICAL_SURFACE"));
        assert!(step.contains("CYLINDRICAL_SURFACE"));
        assert!(step.contains("TOROIDAL_SURFACE"));
        assert!(step.contains("CONICAL_SURFACE"));

        // Full export should carry standalone 1D entities as a secondary wireframe rep.
        assert!(step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"));
    }

    #[test]
    fn hfss_wire_curve_set_references_geometry_not_edge_curve() {
        use std::collections::HashMap;

        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse"); let brep_t = brep; let brep = to_flat(&brep_t);
        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        let mut entity_by_id: HashMap<u64, String> = HashMap::new();
        for line in step.lines() {
            let line = line.trim();
            if !line.starts_with('#') {
                continue;
            }
            let Some((lhs, rhs)) = line.split_once('=') else {
                continue;
            };
            let Some(id_str) = lhs.strip_prefix('#') else {
                continue;
            };
            let Ok(id) = id_str.parse::<u64>() else {
                continue;
            };
            let kind = rhs
                .split_once('(')
                .map(|(k, _)| k.trim().to_string())
                .unwrap_or_default();
            entity_by_id.insert(id, kind);
        }

        let mut curve_set_refs = Vec::new();
        for line in step.lines() {
            let line = line.trim();
            if !line.contains("GEOMETRIC_CURVE_SET") {
                continue;
            }
            let Some(start) = line.find(",(") else {
                continue;
            };
            let Some(end) = line.rfind("))") else {
                continue;
            };
            let refs_text = &line[start + 2..end];
            for token in refs_text.split(',') {
                let tok = token.trim();
                if let Some(id_str) = tok.strip_prefix('#')
                    && let Ok(id) = id_str.parse::<u64>()
                {
                    curve_set_refs.push(id);
                }
            }
        }

        assert!(!curve_set_refs.is_empty(), "expected non-empty GEOMETRIC_CURVE_SET refs");
        for id in curve_set_refs {
            let kind = entity_by_id.get(&id).cloned().unwrap_or_default();
            assert_ne!(
                kind,
                "EDGE_CURVE",
                "GEOMETRIC_CURVE_SET should reference geometric curves, got EDGE_CURVE #{id}"
            );
            assert_ne!(
                kind,
                "SURFACE_CURVE",
                "GEOMETRIC_CURVE_SET should avoid SURFACE_CURVE refs, got SURFACE_CURVE #{id}"
            );
        }
    }

    #[test]
    fn round_trips_sphere_and_cone_surfaces() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse"); let brep_t = brep; let brep = to_flat(&brep_t);

        // Find the original cone half-angle and radius for comparison
        let mut orig_cone_angle = 0.0f64;
        let mut orig_cone_radius = 0.0f64;
        for surface in &brep.geom.surfaces {
            if let Surface3::Cone(c) = surface {
                orig_cone_angle = c.half_angle_rad;
                orig_cone_radius = c.radius;
            }
        }
        assert!(orig_cone_angle > 0.0, "should find a cone in hfss.step");

        // Use non-strict mode for complex geometry round-trip tests
        let step = StepWriter::write_string_with_options(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                ..Default::default()
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("re-exported STEP should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);

        // Count faces with each surface type and verify cone parameters survive round-trip
        let mut sphere_count = 0usize;
        let mut cone_count = 0usize;
        for sid in reparsed.geom.face_surface.iter().flatten() {
            match reparsed.geom.surfaces.get(*sid) {
                Some(Surface3::Sphere(_)) => sphere_count += 1,
                Some(Surface3::Cone(c)) => {
                    cone_count += 1;
                    assert!(
                        (c.half_angle_rad - orig_cone_angle).abs() < 1e-6,
                        "cone half-angle drifted: original={} reparsed={}",
                        orig_cone_angle,
                        c.half_angle_rad,
                    );
                    assert!(
                        (c.radius - orig_cone_radius).abs() < 1e-6,
                        "cone radius drifted: original={} reparsed={}",
                        orig_cone_radius,
                        c.radius,
                    );
                }
                _ => {}
            }
        }
        assert!(
            sphere_count >= 1,
            "expected at least 1 sphere face after round-trip, got {}",
            sphere_count
        );
        assert!(
            cone_count >= 1,
            "expected at least 1 cone face after round-trip, got {}",
            cone_count
        );
    }

    #[test]
    fn creator_style_hfss_keeps_sphere_cylinder_cone_as_solids() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse"); let brep_t = brep; let brep = to_flat(&brep_t);
        let original_solid_count = brep.solids.len();

        // Match creator save path: standard write_string defaults.
        let step = StepWriter::write_string(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        // Solid topology must remain in exported STEP.
        assert!(
            step.contains("MANIFOLD_SOLID_BREP"),
            "export should contain solid entities"
        );

        let reparsed = StepReader::parse_string(&step).expect("re-exported STEP should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);
        assert!(
            !reparsed.solids.is_empty(),
            "re-exported model should contain solids"
        );
        assert_eq!(
            reparsed.solids.len(),
            original_solid_count,
            "solid count should remain stable for creator-style save"
        );

        let mut sphere_faces = 0usize;
        let mut cylinder_faces = 0usize;
        let mut cone_faces = 0usize;
        let mut sphere_solids = 0usize;
        let mut cylinder_solids = 0usize;
        let mut cone_solids = 0usize;
        let mut face_flat_idx = 0usize;
        for solid in &reparsed.solids {
            let mut solid_has_sphere = false;
            let mut solid_has_cylinder = false;
            let mut solid_has_cone = false;
            for shell in &solid.shells {
                for _face in &shell.faces {
                    if let Some(Some(surface_idx)) = reparsed.geom.face_surface.get(face_flat_idx) {
                        match reparsed.geom.surfaces.get(*surface_idx) {
                            Some(Surface3::Sphere(_)) => solid_has_sphere = true,
                            Some(Surface3::Cylinder(_)) => solid_has_cylinder = true,
                            Some(Surface3::Cone(_)) => solid_has_cone = true,
                            _ => {}
                        }
                    }
                    face_flat_idx += 1;
                }
            }
            if solid_has_sphere {
                sphere_solids += 1;
            }
            if solid_has_cylinder {
                cylinder_solids += 1;
            }
            if solid_has_cone {
                cone_solids += 1;
            }
        }
        for sid in reparsed.geom.face_surface.iter().flatten() {
            match reparsed.geom.surfaces.get(*sid) {
                Some(Surface3::Sphere(_)) => sphere_faces += 1,
                Some(Surface3::Cylinder(_)) => cylinder_faces += 1,
                Some(Surface3::Cone(_)) => cone_faces += 1,
                _ => {}
            }
        }
        assert!(sphere_faces > 0, "sphere should remain analytic after save");
        assert!(cylinder_faces > 0, "cylinder should remain analytic after save");
        assert!(cone_faces > 0, "cone should remain analytic after save");
        assert!(sphere_solids > 0, "sphere should belong to at least one solid");
        assert!(
            cylinder_solids > 0,
            "cylinder should belong to at least one solid"
        );
        assert!(cone_solids > 0, "cone should belong to at least one solid");
    }

    #[test]
    fn exports_ellipsoid_surface_emits_semantic_tag() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse"); let brep_t = brep; let mut brep = to_flat(&brep_t);
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Ellipsoid(rcad_kernel::EllipsoidalSurface {
            center: DVec3::new(0.5, 0.5, 0.5),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 2.0,
            radius_y: 1.5,
            radius_z: 1.0,
        });

        let step = StepWriter::write_string_with_options(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                ..Default::default()
            },
        );
        assert!(step.contains("B_SPLINE_SURFACE_WITH_KNOTS"));
        assert!(step.contains("RCAD_ELLIPSOID"));

        let reparsed = StepReader::parse_string(&step).expect("ellipsoid fallback STEP should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);
        let ellipsoid_surfaces = reparsed
            .geom
            .surfaces
            .iter()
            .filter(|surface| matches!(surface, Surface3::Ellipsoid(_)))
            .count();
        assert!(
            ellipsoid_surfaces > 0,
            "expected at least one reparsed ellipsoid surface"
        );
    }

    #[test]
    fn exports_helicoid_surface_emits_semantic_tag() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse"); let brep_t = brep; let mut brep = to_flat(&brep_t);
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Helicoid(rcad_kernel::HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 3.0,
        });

        let step = StepWriter::write_string_with_options(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                ..Default::default()
            },
        );
        assert!(step.contains("B_SPLINE_SURFACE_WITH_KNOTS"));
        assert!(step.contains("RCAD_HELICOID"));

        let reparsed = StepReader::parse_string(&step).expect("helicoid fallback STEP should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);
        let helicoid_surfaces = reparsed
            .geom
            .surfaces
            .iter()
            .filter(|surface| matches!(surface, Surface3::Helicoid(_)))
            .count();
        assert!(
            helicoid_surfaces > 0,
            "expected at least one reparsed helicoid surface"
        );
    }

    #[test]
    fn exports_coons_surface_via_bspline_fallback() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse"); let brep_t = brep; let mut brep = to_flat(&brep_t);
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Coons(rcad_kernel::CoonsSurface {
            south: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::X,
            })),
            north: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(0.0, 1.0, 1.0),
                direction: DVec3::X,
            })),
            west: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
            east: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(1.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
        });

        let step = StepWriter::write_string_with_options(
            &brep.to_topods(),
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                ..Default::default()
            },
        );
        assert!(step.contains("B_SPLINE_SURFACE_WITH_KNOTS"));
        assert!(!step.contains("RCAD_COONS"));

        let reparsed = StepReader::parse_string(&step).expect("Coons fallback STEP should parse"); let reparsed_t = reparsed; let reparsed = to_flat(&reparsed_t);
        let bspline_surfaces = reparsed
            .geom
            .surfaces
            .iter()
            .filter(|surface| matches!(surface, Surface3::BSpline(_)))
            .count();
        assert!(
            bspline_surfaces > 0,
            "expected at least one reparsed bspline surface from Coons fallback"
        );
