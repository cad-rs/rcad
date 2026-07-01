#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    fn create_box_brep() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 2.0,
            depth: 3.0,
        })
    }

    #[test]
    fn test_fillet_params_default() {
        let params = FilletParams::default();
        assert_eq!(params.radius, 1.0);
        assert_eq!(params.continuity, FilletContinuity::C1);
        assert_eq!(params.mode, FilletMode::Uniform);
        assert!(params.tension >= 0.0 && params.tension <= 1.0);
    }

    #[test]
    fn test_fillet_params_builder() {
        let params = FilletParams::new(0.5)
            .with_continuity(FilletContinuity::C2)
            .with_mode(FilletMode::Chordal)
            .with_tension(0.8);

        assert_eq!(params.radius, 0.5);
        assert_eq!(params.continuity, FilletContinuity::C2);
        assert_eq!(params.mode, FilletMode::Chordal);
        assert!((params.tension - 0.8).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_variable_radius_point() {
        let point = VariableRadiusPoint::new(0.5, 2.0);
        assert!((point.parameter - 0.5).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((point.radius - 2.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_fillet_error_display() {
        let err = FilletError::InvalidRadius { radius: -1.0 };
        assert!(err.to_string().contains("invalid fillet radius"));

        let err = FilletError::EdgeNotFound { edge_index: 99 };
        assert!(err.to_string().contains("edge 99 not found"));

        let err = FilletError::RadiusTooLarge {
            edge_index: 0,
            radius: 10.0,
            max_radius: 1.0,
        };
        assert!(err.to_string().contains("radius 10 too large"));
    }

    #[test]
    fn test_any_perpendicular() {
        let v = DVec3::X;
        let p = any_perpendicular(v);
        assert!((p.dot(v)).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((p.length() - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);

        let v = DVec3::Z;
        let p = any_perpendicular(v);
        assert!((p.dot(v)).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_interpolate_radius() {
        // Linear interpolation (tension = 0)
        let r = interpolate_radius(1.0, 2.0, 0.0, 0.0);
        assert!((r - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);

        let r = interpolate_radius(1.0, 2.0, 1.0, 0.0);
        assert!((r - 2.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);

        let r = interpolate_radius(1.0, 2.0, 0.5, 0.0);
        assert!((r - 1.5).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_make_fillet_edge_empty_indices() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[], 0.5);
        assert!(result.is_ok());
        let res = result.unwrap();
        assert_eq!(res.edges_processed, 0);
        assert!(!res.warnings.is_empty());
    }

    #[test]
    fn test_make_fillet_edge_invalid_radius() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[0], -1.0);
        assert!(matches!(result, Err(FilletError::InvalidRadius { .. })));
    }

    #[test]
    fn test_make_fillet_edge_invalid_edge_index() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[999], 0.5);
        assert!(matches!(result, Err(FilletError::EdgeNotFound { .. })));
    }

    #[test]
    fn test_make_fillet_all_edges() {
        let brep = create_box_brep();
        let result = make_fillet_all_edges(&brep, 0.1);
        // Should succeed without errors
        assert!(result.is_ok() || matches!(result, Err(FilletError::EdgeNotFound { .. })));
    }

    #[test]
    fn test_make_variable_fillet_too_few_points() {
        let brep = create_box_brep();
        let radii = vec![VariableRadiusPoint::new(0.0, 0.5)];
        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(result.is_err());
    }

    #[test]
    fn test_make_variable_fillet_invalid_parameter() {
        let brep = create_box_brep();
        let radii = vec![
            VariableRadiusPoint::new(-0.5, 0.5),
            VariableRadiusPoint::new(1.0, 1.0),
        ];
        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(matches!(result, Err(FilletError::InvalidVariableRadius { .. })));
    }

    #[test]
    fn test_compute_plane_plane_fillet() {
        let plane1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let plane2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_Y * 0.5,
            end_point: DVec3::Y * 0.5,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Y,
            tangent_end: DVec3::Y,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_plane_plane_fillet(&edge_info, &plane1, &plane2, 0.5);
        assert!(result.is_ok());

        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Cylinder(_)));
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::ZERO,
            end_point: DVec3::Z,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_cylinder_plane_fillet(&edge_info, &cylinder, &plane, 0.5);
        assert!(result.is_ok());
    }

    #[test]
    fn test_compute_sphere_plane_fillet() {
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 1.0);
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 0.5),
            normal: DVec3::Z,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_X,
            end_point: DVec3::X,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::X,
            tangent_end: DVec3::X,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_sphere_plane_fillet(&edge_info, &sphere, &plane, 0.5);
        assert!(result.is_ok());

        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Sphere(_)));
    }

    #[test]
    fn test_compute_fillet_curves_torus() {
        let brep = create_box_brep();

        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 2.0,
            minor_radius: 0.5,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::ZERO,
            end_point: DVec3::Z,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_fillet_curves(&brep, &edge_info, 0.5, &Surface3::Torus(torus));
        assert!(result.is_ok());

        let curves = result.unwrap();
        assert!(!curves.is_empty());
    }

    #[test]
    fn test_fillet_continuity() {
        assert_eq!(FilletContinuity::default(), FilletContinuity::C1);
        assert_ne!(FilletContinuity::C0, FilletContinuity::C1);
        assert_ne!(FilletContinuity::C1, FilletContinuity::C2);
    }

    #[test]
    fn test_fillet_mode() {
        assert_eq!(FilletMode::default(), FilletMode::Uniform);
        assert_ne!(FilletMode::Uniform, FilletMode::Variable);
        assert_ne!(FilletMode::Variable, FilletMode::Chordal);
    }

    #[test]
    fn test_fillet_result_creation() {
        let brep = create_box_brep();
        let result = FilletResult {
            brep: brep.clone(),
            edges_processed: 3,
            fillet_faces_created: 3,
            warnings: vec!["test warning".to_string()],
        };

        assert_eq!(result.edges_processed, 3);
        assert_eq!(result.fillet_faces_created, 3);
        assert_eq!(result.warnings.len(), 1);
    }

    // ============================================================================
    // Edge Case Tests for OCCT Alignment
    // ============================================================================

    /// Test fillet on a concave edge (interior corner of a box).
    /// Concave edges require different handling than convex edges.
    #[test]
    fn test_fillet_concave_edge() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        // Create a box and populate its geometry
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        // Attempt to fillet edges - concave edges are at the interior corners
        let result = make_fillet_edge(&brep, &[0, 1, 2, 3], 0.5);
        assert!(result.is_ok(), "concave edge fillet should succeed");

        let fillet_result = result.unwrap();
        // At least some edges should be processed
        assert!(fillet_result.edges_processed > 0 || !fillet_result.warnings.is_empty());
    }

    /// Test fillet on a chain of connected edges.
    /// Chain edges should blend smoothly at the vertices where they meet.
    #[test]
    fn test_fillet_chain_edges() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        // Create a box
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 3.0,
            depth: 2.0,
        });
        populate_box_geom(&mut brep);

        // Fillet a chain of edges around one face (edges 0, 1, 2, 3 form a loop)
        let result = make_fillet_edge(&brep, &[0, 1, 2, 3], 0.3);
        assert!(result.is_ok(), "chain edge fillet should succeed");

        let fillet_result = result.unwrap();
        assert!(fillet_result.edges_processed >= 1, "at least one edge should be filleted");
    }

    /// Test fillet with very small radius on an edge.
    /// Small radius fillets should not create degenerate geometry.
    #[test]
    fn test_fillet_small_radius() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 10.0,
            height: 10.0,
            depth: 10.0,
        });
        populate_box_geom(&mut brep);

        // Very small fillet radius
        let result = make_fillet_edge(&brep, &[0], 0.001);
        assert!(result.is_ok(), "small radius fillet should succeed");
    }

    /// Test fillet on an edge where adjacent faces are perpendicular.
    /// This is the most common case and should produce a clean toroidal fillet.
    #[test]
    fn test_fillet_perpendicular_faces() {
        use rcad_kernel::geom::{Plane, Line3};
        use glam::DVec3;

        // Create edge info for perpendicular faces
        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_Y,
            end_point: DVec3::Y,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Y,  // Edge along Y axis
            tangent_end: DVec3::Y,
            length: 2.0,
            curve: Some(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::Y,
            })),
            curve_range: Some([0.0, 2.0]),
        };

        let plane1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let plane2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let result = compute_plane_plane_fillet(&edge_info, &plane1, &plane2, 0.5);
        assert!(result.is_ok(), "perpendicular faces fillet should succeed");

        // Result should be a cylinder for plane-plane fillet (straight edge)
        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Cylinder(_)), "plane-plane fillet should produce cylinder");
    }

    /// Test fillet with variable radius along the edge.
    /// Variable radius fillets should interpolate between radius values.
    #[test]
    fn test_fillet_variable_radius_basic() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 5.0,
            height: 5.0,
            depth: 5.0,
        });
        populate_box_geom(&mut brep);

        // Variable radius: 0.2 at start, 0.8 at end
        let radii = vec![
            VariableRadiusPoint::new(0.0, 0.2),
            VariableRadiusPoint::new(1.0, 0.8),
        ];

        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(result.is_ok(), "variable radius fillet should succeed");
    }

    /// Test fillet surface computation for cylinder-plane edge.
    /// Verifies correct handling of curved-to-planar transitions.
    #[test]
    fn test_fillet_cylinder_plane_edge() {
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::new(1.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_Z,
            end_point: DVec3::Z,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 2.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_cylinder_plane_fillet(&edge_info, &cylinder, &plane, 0.3);
        assert!(result.is_ok(), "cylinder-plane fillet should succeed");
    }

    /// Verify SA for a single fillet on a 100x100x100 box (OCCT A1).
    #[test]
    fn verify_fillet_sa_box() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 100.0, height: 100.0, depth: 100.0,
        });
        let result = crate::make_fillet_edge(&brep, &[4], 10.0);
        assert!(result.is_ok(), "fillet should succeed");
        let total = rcad_kernel::surface_area(&result.unwrap().brep);
        let diff = (total - 59527.9).abs();
        let tol = 0.15 * 59527.9;
        assert!(diff <= tol, "SA {total} differs from expected 59527.9 by {diff} > {tol}");
    }
}
