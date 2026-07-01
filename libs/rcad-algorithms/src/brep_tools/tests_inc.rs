mod tests {
    use super::*;
    use crate::tolerance::*;
    use rcad_kernel::PrimitiveSolid;
    use std::f64::consts::PI;

    fn make_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 2.0,
            depth: 3.0,
        })
    }

    // 鈹€鈹€ I/O Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_write_brep_to_string() {
        let brep = make_box();
        let json = write_brep_to_string(&brep).unwrap();
        assert!(json.contains("vertices"));
        assert!(json.contains("edges"));
        assert!(json.contains("solids"));
    }

    #[test]
    fn test_read_brep_from_string() {
        let brep = make_box();
        let json = write_brep_to_string(&brep).unwrap();
        let restored = read_brep_from_string(&json).unwrap();

        assert_eq!(brep.vertices.len(), restored.vertices.len());
        assert_eq!(brep.edges.len(), restored.edges.len());
        assert_eq!(brep.solids.len(), restored.solids.len());
    }

    #[test]
    fn test_read_invalid_json() {
        let result = read_brep_from_string("not valid json");
        assert!(result.is_err());
    }

    #[test]
    fn test_write_and_read_roundtrip() {
        let original = make_box();
        let json = write_brep_to_string(&original).unwrap();
        let restored = read_brep_from_string(&json).unwrap();

        // Check vertices match
        for (orig, rest) in original.vertices.iter().zip(restored.vertices.iter()) {
            assert!((orig.point - rest.point).length() < TOLERANCE_COORD_SUB);
        }
    }

    // 鈹€鈹€ Transformation Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_transform_shape_translation() {
        let mut brep = make_box();
        let original_vertex = brep.vertices[0].point;

        let translation = DAffine3::from_translation(DVec3::new(5.0, 0.0, 0.0));
        transform_shape(&mut brep, translation);

        let expected = original_vertex + DVec3::new(5.0, 0.0, 0.0);
        assert!((brep.vertices[0].point - expected).length() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_transform_shape_rotation() {
        let mut brep = make_box();
        // Use vertex 1 which is at (width, 0, 0), not at origin
        let original_vertex = brep.vertices[1].point;

        // Rotate 90 degrees around Z axis
        let rotation = DAffine3::from_axis_angle(DVec3::Z, PI / 2.0);
        transform_shape(&mut brep, rotation);

        // The vertex should have moved (vertex 1 is at (1, 0, 0), after rotation it's at (0, 1, 0))
        assert!((brep.vertices[1].point - original_vertex).length() > 0.1);
    }

    #[test]
    fn test_mirror_shape() {
        let mut brep = make_box();
        let original_x = brep.vertices.iter()
            .map(|v| v.point.x)
            .fold(f64::INFINITY, |a, b| a.min(b));

        // Mirror across the YZ plane at x=0
        mirror_shape(&mut brep, DVec3::ZERO, DVec3::X);

        // The minimum X should now be negative (mirrored)
        let new_min_x = brep.vertices.iter()
            .map(|v| v.point.x)
            .fold(f64::INFINITY, |a, b| a.min(b));

        assert!(new_min_x < 0.0);
    }

    #[test]
    fn test_scale_shape() {
        let mut brep = make_box();

        let original_volume = rcad_kernel::volume(&brep);

        // Scale by 2x about the origin
        scale_shape(&mut brep, 2.0, DVec3::ZERO);

        let new_volume = rcad_kernel::volume(&brep);

        // Volume should scale by 2^3 = 8
        assert!((new_volume / original_volume - 8.0).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_rotate_shape() {
        let mut brep = make_box();
        let original_bb = bounding_box(&brep).unwrap();

        // Rotate 90 degrees around Z axis through origin
        rotate_shape(&mut brep, DVec3::ZERO, DVec3::Z, PI / 2.0);

        let new_bb = bounding_box(&brep).unwrap();

        // After 90-degree rotation, the bounding box dimensions should swap
        // Original: dx=1, dy=2 -> after rotation: dx=2, dy=1
        let original_size = original_bb[1] - original_bb[0];
        let new_size = new_bb[1] - new_bb[0];

        // X and Y dimensions should have swapped
        assert!((original_size.x - new_size.y).abs() < TOLERANCE_COORD_SUB);
        assert!((original_size.y - new_size.x).abs() < TOLERANCE_COORD_SUB);
        assert!((original_size.z - new_size.z).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_rotate_about_arbitrary_axis() {
        let mut brep = make_box();

        // Rotate 180 degrees about an axis through the center
        let center = DVec3::new(0.5, 1.0, 1.5);
        rotate_shape(&mut brep, center, DVec3::Z, PI);

        // The box should still have the same volume
        let volume = rcad_kernel::volume(&brep);
        assert!((volume - 6.0).abs() < TOLERANCE_MESH_LEGACY); // 1 * 2 * 3
    }

    // 鈹€鈹€ Shape Type Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_get_shape_type_solid() {
        let brep = make_box();
        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
    }

    #[test]
    fn test_get_shape_type_empty() {
        let brep = BRep::new();
        assert_eq!(get_shape_type(&brep), ShapeType::Empty);
    }

    #[test]
    fn test_get_shape_type_compound() {
        let mut compound = rcad_kernel::topology::Compound::new();
        compound.add_solid(None, rcad_kernel::topology::Solid { shells: vec![] });
        let brep = BRep::from_compound(compound);
        assert_eq!(get_shape_type(&brep), ShapeType::Compound);
    }

    // 鈹€鈹€ Closure Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_is_closed_box() {
        let brep = make_box();
        assert!(is_closed(&brep));
    }

    #[test]
    fn test_is_closed_empty() {
        let brep = BRep::new();
        assert!(!is_closed(&brep));
    }

    #[test]
    fn test_is_closed_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        assert!(is_closed(&brep));
    }

    #[test]
    fn test_is_closed_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        assert!(is_closed(&brep));
    }

    // 鈹€鈹€ Count Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_count_faces() {
        let brep = make_box();
        assert_eq!(count_faces(&brep), 6); // Box has 6 faces
    }

    #[test]
    fn test_count_edges() {
        let brep = make_box();
        assert_eq!(count_edges(&brep), 12); // Box has 12 edges
    }

    #[test]
    fn test_count_vertices() {
        let brep = make_box();
        assert_eq!(count_vertices(&brep), 8); // Box has 8 vertices
    }

    #[test]
    fn test_count_shells() {
        let brep = make_box();
        assert_eq!(count_shells(&brep), 1); // Box has 1 shell
    }

    // 鈹€鈹€ Bounding Box Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_bounding_box_box() {
        let brep = make_box();
        let bb = bounding_box(&brep).unwrap();

        assert!((bb[0].x - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[0].y - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[0].z - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].x - 1.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].y - 2.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].z - 3.0).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_bounding_box_empty() {
        let brep = BRep::new();
        assert!(bounding_box(&brep).is_none());
    }

    #[test]
    fn test_bounding_box_sphere() {
        // Note: Sphere primitive only has 2 vertices (poles at y=+r and y=-r)
        // The bounding box based on vertices will only cover the poles
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let bb = bounding_box(&brep).unwrap();

        // The sphere has vertices at (0, +r, 0) and (0, -r, 0)
        // So bounding box is: min=(0, -1, 0), max=(0, 1, 0)
        assert!((bb[0].y - (-1.0)).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].y - 1.0).abs() < TOLERANCE_COORD_SUB);
        // X and Z are 0 because only pole vertices exist
        assert!((bb[0].x - 0.0).abs() < TOLERANCE_COORD_SUB);
        assert!((bb[1].x - 0.0).abs() < TOLERANCE_COORD_SUB);
    }

    // 鈹€鈹€ Wire Query Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_get_outer_wire() {
        let brep = make_box();
        let wire = get_outer_wire(&brep, 0).unwrap();
        assert_eq!(wire.edges.len(), 4); // Each box face is a quad
    }

    #[test]
    fn test_get_outer_wire_invalid_index() {
        let brep = make_box();
        let result = get_outer_wire(&brep, 100);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_inner_wires_empty() {
        let brep = make_box();
        let inner = get_inner_wires(&brep, 0).unwrap();
        assert!(inner.is_empty()); // Box faces have no holes
    }

    // 鈹€鈹€ Shape Type Display Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_shape_type_display() {
        assert_eq!(format!("{}", ShapeType::Solid), "Solid");
        assert_eq!(format!("{}", ShapeType::Face), "Face");
        assert_eq!(format!("{}", ShapeType::Compound), "Compound");
        assert_eq!(format!("{}", ShapeType::Empty), "Empty");
    }

    // 鈹€鈹€ Error Display Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_error_display() {
        let err = BRepToolsError::InvalidIndex {
            kind: "face",
            index: 10,
            max: 5,
        };
        assert!(format!("{}", err).contains("Invalid face index"));

        let err = BRepToolsError::MissingGeometry {
            kind: "surface",
            index: 5,
        };
        assert!(format!("{}", err).contains("Missing surface"));
    }

    // 鈹€鈹€ Integration Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn test_transform_and_serialize() {
        let mut brep = make_box();

        // Apply transformation
        scale_shape(&mut brep, 2.0, DVec3::ZERO);

        // Serialize
        let json = write_brep_to_string(&brep).unwrap();

        // Deserialize and verify
        let restored = read_brep_from_string(&json).unwrap();
        let restored_volume = rcad_kernel::volume(&restored);

        // Volume should be 6.0 * 8 = 48.0
        assert!((restored_volume - 48.0).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_multiple_transformations() {
        let mut brep = make_box();
        let original_volume = rcad_kernel::volume(&brep);

        // Apply multiple transformations
        rotate_shape(&mut brep, DVec3::ZERO, DVec3::Z, PI / 4.0);
        scale_shape(&mut brep, 1.5, DVec3::ZERO);
        let translation = DAffine3::from_translation(DVec3::new(10.0, 0.0, 0.0));
        transform_shape(&mut brep, translation);

        // Volume should be scaled by 1.5^3 = 3.375
        let new_volume = rcad_kernel::volume(&brep);
        assert!((new_volume / original_volume - 3.375).abs() < TOLERANCE_MESH_LEGACY);

        // Bounding box should be shifted
        let bb = bounding_box(&brep).unwrap();
        assert!(bb[0].x > 5.0); // Should be shifted in positive X
    }

    #[test]
    fn test_sphere_operations() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Check that vertices scale correctly
        // The sphere has vertices at poles: (0, r, 0) and (0, -r, 0)
        let original_y = brep.vertices[0].point.y;

        // Scale sphere by 2x
        scale_shape(&mut brep, 2.0, DVec3::ZERO);

        // Vertex should be scaled by 2
        assert!((brep.vertices[0].point.y - original_y * 2.0).abs() < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn test_cylinder_operations() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Cylinder has 3 faces (top, bottom, side)
        assert_eq!(count_faces(&brep), 3);
    }

    #[test]
    fn test_cone_operations() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 1.0,
            height: 2.0,
        });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Cone has 2 faces (base, side)
        assert_eq!(count_faces(&brep), 2);
    }

    #[test]
    fn test_torus_operations() {
        let brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        assert_eq!(get_shape_type(&brep), ShapeType::Solid);
        assert!(is_closed(&brep));

        // Torus has 1 face
        assert_eq!(count_faces(&brep), 1);
    }
}
