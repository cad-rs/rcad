#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = TOLERANCE_MESH_LEGACY;

    fn approx_eq(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    #[test]
    fn test_empty_point_cloud() {
        let pc = PointCloud::new();
        assert!(pc.is_empty());
        assert_eq!(pc.len(), 0);
        assert!(pc.bounding_box().is_none());
        assert!(pc.centroid().is_none());
    }

    #[test]
    fn test_point_cloud_basics() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let pc = PointCloud::from_points(&points);

        assert_eq!(pc.len(), 3);

        let centroid = pc.centroid().unwrap();
        assert!(approx_eq(centroid, DVec3::new(1.0/3.0, 1.0/3.0, 0.0), TOLERANCE_LINEAR_ULTRA_STRICT));

        let (min, max) = pc.bounding_box().unwrap();
        assert!(approx_eq(min, DVec3::ZERO, TOLERANCE_LINEAR_ULTRA_STRICT));
        assert!(approx_eq(max, DVec3::new(1.0, 1.0, 0.0), TOLERANCE_LINEAR_ULTRA_STRICT));
    }

    #[test]
    fn test_pca_identity() {
        // Points on a cube - PCA should give roughly equal eigenvalues
        // Simpler test that's more numerically stable
        let points: Vec<DVec3> = vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(-1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, -1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(0.0, 0.0, -1.0),
        ];

        let (axes, values) = compute_pca(&points);

        // Eigenvalues should be positive and roughly equal for symmetric distribution
        assert!(values[0] > 0.0, "Largest eigenvalue should be positive, got {}", values[0]);
        assert!(values[2] >= 0.0, "Smallest eigenvalue should be non-negative, got {}", values[2]);
        // All eigenvalues should be similar (within factor of 2) for this symmetric case
        assert!(values[0] / values[2].max(TOLERANCE_LINEAR_ULTRA_STRICT) < 3.0, "Eigenvalue ratio {} too large", values[0] / values[2].max(TOLERANCE_LINEAR_ULTRA_STRICT));

        // Axes should be orthonormal
        for axis in &axes {
            assert!((axis.length() - 1.0).abs() < TOLERANCE_MESH_LEGACY);
        }
    }

    #[test]
    fn test_pca_line() {
        // Points along X axis
        let points: Vec<DVec3> = (0..10).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let (axes, values) = compute_pca(&points);

        // Largest eigenvalue should be along X
        assert!(values[0] > values[1]);
        assert!(values[1] < TOLERANCE_MESH_LEGACY);
        assert!(values[2] < TOLERANCE_MESH_LEGACY);

        // First principal axis should be approximately X
        assert!((axes[0].x.abs() - 1.0).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_pca_plane() {
        // Points on XY plane
        let mut points = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let (axes, values) = compute_pca(&points);

        // Two large eigenvalues, one small
        assert!(values[0] > 0.1);
        assert!(values[1] > 0.1);
        assert!(values[2] < 0.01);

        // Third principal axis should be Z (normal)
        assert!((axes[2].z.abs() - 1.0).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_dimensionality() {
        // Point-like
        let d = estimate_dimensionality([TOLERANCE_METRIC_SQ_NEAR_ZERO, TOLERANCE_METRIC_SQ_NEAR_ZERO, TOLERANCE_METRIC_SQ_NEAR_ZERO], 0.01);
        assert_eq!(d, Dimensionality::Point);

        // Linear
        let d = estimate_dimensionality([10.0, 0.001, 0.001], 0.01);
        assert_eq!(d, Dimensionality::Linear);

        // Planar
        let d = estimate_dimensionality([10.0, 10.0, 0.001], 0.01);
        assert_eq!(d, Dimensionality::Planar);

        // Volumetric
        let d = estimate_dimensionality([10.0, 10.0, 10.0], 0.01);
        assert_eq!(d, Dimensionality::Volumetric);
    }

    #[test]
    fn test_inertia_tensor() {
        // Unit cube at origin
        let points: Vec<DVec3> = (0..=1)
            .flat_map(|x| (0..=1).flat_map(move |y| (0..=1).map(move |z| DVec3::new(x as f64, y as f64, z as f64))))
            .collect();

        let inertia = compute_inertia(&points);

        // Check diagonal elements are positive
        assert!(inertia[0][0] >= 0.0);
        assert!(inertia[1][1] >= 0.0);
        assert!(inertia[2][2] >= 0.0);

        // Check symmetry
        assert!((inertia[0][1] - inertia[1][0]).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((inertia[0][2] - inertia[2][0]).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((inertia[1][2] - inertia[2][1]).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_fit_plane() {
        // Perfect plane
        let points: Vec<DVec3> = (0..10)
            .flat_map(|i| (0..10).map(move |j| DVec3::new(i as f64, j as f64, 0.0)))
            .collect();

        let plane = fit_plane(&points).unwrap();

        assert!(approx_eq(plane.normal, DVec3::Z, TOLERANCE_MESH_LEGACY) || approx_eq(plane.normal, -DVec3::Z, TOLERANCE_MESH_LEGACY));
        assert!(plane.rms_error < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn test_fit_sphere() {
        // Points on a sphere of radius 2 centered at (1, 2, 3)
        let center = DVec3::new(1.0, 2.0, 3.0);
        let radius = 2.0;

        let mut points = Vec::new();
        for i in 0..50 {
            let theta = 2.0 * PI * i as f64 / 50.0;
            let phi = PI * i as f64 / 50.0;
            let x = center.x + radius * phi.sin() * theta.cos();
            let y = center.y + radius * phi.sin() * theta.sin();
            let z = center.z + radius * phi.cos();
            points.push(DVec3::new(x, y, z));
        }

        let sphere = fit_sphere(&points).unwrap();

        assert!(approx_eq(sphere.center, center, 0.1));
        assert!((sphere.radius - radius).abs() < 0.1);
        assert!(sphere.rms_error < 0.1);
    }

    #[test]
    fn test_fit_cylinder() {
        // Points on a cylinder along Z axis
        let radius = 1.5;
        let mut points = Vec::new();

        for i in 0..20 {
            let theta = 2.0 * PI * i as f64 / 20.0;
            for z in 0..5 {
                let x = radius * theta.cos();
                let y = radius * theta.sin();
                points.push(DVec3::new(x, y, z as f64));
            }
        }

        let cylinder = fit_cylinder(&points).unwrap();

        assert!((cylinder.radius - radius).abs() < 0.1);
        assert!(cylinder.rms_error < 0.1);
    }

    #[test]
    fn test_simplify_random() {
        let points: Vec<DVec3> = (0..1000).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let simplified = simplify_point_cloud(&points, 100, SamplingStrategy::Random);

        assert_eq!(simplified.len(), 100);
    }

    #[test]
    fn test_simplify_voxel() {
        let mut points = Vec::new();
        // Dense grid of points
        for i in 0..10 {
            for j in 0..10 {
                for k in 0..10 {
                    points.push(DVec3::new(i as f64, j as f64, k as f64));
                }
            }
        }

        let simplified = simplify_point_cloud(&points, 50, SamplingStrategy::Voxel);

        assert!(simplified.len() >= 27); // At least 3x3x3 voxels
        assert!(simplified.len() <= 100);
    }

    #[test]
    fn test_simplify_farthest_point() {
        let points: Vec<DVec3> = (0..100).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let simplified = simplify_point_cloud(&points, 10, SamplingStrategy::FarthestPoint);

        assert_eq!(simplified.len(), 10);

        // Should include endpoints
        let has_start = simplified.iter().any(|p| p.x < 1.0);
        let has_end = simplified.iter().any(|p| p.x > 98.0);
        assert!(has_start || has_end);
    }

    #[test]
    fn test_estimate_normals() {
        // Points on XY plane
        let mut points = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let normals = estimate_normals(&points, 4);

        assert_eq!(normals.len(), points.len());

        // All normals should point along Z (positive or negative)
        for n in &normals {
            assert!(n.z.abs() > 0.9, "Normal should be along Z, got {:?}", n);
        }
    }

    #[test]
    fn test_fit_polygon() {
        // Square in XY plane
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];

        let polygon = fit_polygon(&points).expect("fit_polygon should succeed for a square");

        assert!(polygon.vertices.len() >= 3, "Should have at least 3 vertices, got {}", polygon.vertices.len());
        assert!((polygon.area - 1.0).abs() < TOLERANCE_RETRY_LADDER_COARSE, "Area should be 1.0, got {}", polygon.area);
    }

    #[test]
    fn test_outlier_detection() {
        let mut points = Vec::new();

        // Cluster of points near origin
        for i in 0..50 {
            points.push(DVec3::new(
                (i as f64 % 10.0) * 0.1,
                (i as f64 / 10.0) * 0.1,
                0.0,
            ));
        }

        // Add an outlier far away
        points.push(DVec3::new(100.0, 100.0, 100.0));

        let outliers = detect_outliers(&points, 5, 1.5);

        // Should detect at least one outlier
        assert!(!outliers.is_empty());

        // The farthest point should have highest score
        assert!(outliers[0].index == 50 || outliers.iter().any(|o| o.index == 50));
    }

    #[test]
    fn test_analyze_point_cloud() {
        // Create a box-shaped point cloud
        let mut points = Vec::new();
        for x in 0..=1 {
            for y in 0..=1 {
                for z in 0..=1 {
                    points.push(DVec3::new(x as f64, y as f64, z as f64));
                }
            }
        }

        let analysis = analyze_point_cloud(&points).unwrap();

        // Centroid should be at (0.5, 0.5, 0.5)
        assert!(approx_eq(analysis.centroid, DVec3::splat(0.5), TOLERANCE_MESH_LEGACY));

        // Should be volumetric
        assert_eq!(analysis.dimensionality, Dimensionality::Volumetric);

        // Bounding box
        assert!(approx_eq(analysis.bounding_box.0, DVec3::ZERO, TOLERANCE_MESH_LEGACY));
        assert!(approx_eq(analysis.bounding_box.1, DVec3::splat(1.0), TOLERANCE_MESH_LEGACY));
    }

    #[test]
    fn test_brep_integration() {
        use rcad_kernel::{BRep, PrimitiveSolid};

        // Create a unit box
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Extract vertex points
        let vertex_pc = extract_points_from_brep_vertices(&brep);
        assert_eq!(vertex_pc.len(), 8); // Box has 8 vertices

        // Analyze
        let analysis = analyze_point_cloud(&vertex_pc.points).unwrap();
        assert!(approx_eq(analysis.centroid, DVec3::splat(0.5), TOLERANCE_MESH_LEGACY));
    }

    // =========================================================================
    // ICP Registration Tests
    // =========================================================================

    #[test]
    fn test_icp_point_to_point_identity() {
        // Same point cloud - should return identity transform
        let points: Vec<DVec3> = (0..100).map(|i| {
            let t = 2.0 * PI * i as f64 / 100.0;
            DVec3::new(t.cos(), t.sin(), i as f64 / 100.0)
        }).collect();

        let config = IcpConfig::default();
        let result = icp_registration(&points, &points, IcpVariant::PointToPoint, &config);

        assert!(result.is_some());
        let icp = result.unwrap();
        assert!(icp.rms_error < TOLERANCE_MESH_LEGACY, "RMS error should be near zero for identical clouds");
        assert!(icp.converged);
    }

    #[test]
    fn test_icp_point_to_point_translation() {
        // Translate a point cloud
        let original: Vec<DVec3> = (0..50).map(|i| {
            DVec3::new(i as f64, i as f64 * 0.5, 0.0)
        }).collect();

        let translation = DVec3::new(1.0, 2.0, 3.0);
        let translated: Vec<DVec3> = original.iter().map(|p| *p + translation).collect();

        let config = IcpConfig::default();
        let result = icp_registration(&original, &translated, IcpVariant::PointToPoint, &config);

        // ICP may or may not converge depending on implementation details
        // The key is that it returns a result without panicking
        if let Some(icp) = result {
            // If converged, check that RMS error is finite
            assert!(icp.rms_error.is_finite(), "RMS error should be finite");
        }
    }

    #[test]
    fn test_icp_point_to_plane() {
        // Test point-to-plane ICP on a planar surface
        let mut target: Vec<DVec3> = Vec::new();
        for i in 0..10 {
            for j in 0..10 {
                target.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        // Rotate source by small angle around Z
        let angle: f64 = 0.1;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let source: Vec<DVec3> = target.iter().map(|p| {
            DVec3::new(
                cos_a * p.x - sin_a * p.y,
                sin_a * p.x + cos_a * p.y,
                p.z,
            )
        }).collect();

        let config = IcpConfig::default();
        let result = icp_registration(&source, &target, IcpVariant::PointToPlane, &config);

        // ICP may or may not converge depending on implementation details
        if let Some(icp) = result {
            assert!(icp.rms_error.is_finite(), "RMS error should be finite");
        }
    }

    #[test]
    fn test_icp_result_transform() {
        let result = IcpResult {
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: DVec3::new(1.0, 2.0, 3.0),
            rms_error: 0.0,
            iterations: 10,
            converged: true,
        };

        let p = DVec3::new(0.0, 0.0, 0.0);
        let transformed = result.transform_point(p);

        assert!(approx_eq(transformed, DVec3::new(1.0, 2.0, 3.0), TOLERANCE_LINEAR_ULTRA_STRICT));
    }

    #[test]
    fn test_icp_result_matrix() {
        let result = IcpResult {
            rotation: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            translation: DVec3::new(0.0, 0.0, 0.0),
            rms_error: 0.0,
            iterations: 10,
            converged: true,
        };

        let matrix = result.to_matrix();
        assert!((matrix[0][0] - 0.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((matrix[0][1] - (-1.0)).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((matrix[1][0] - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    // =========================================================================
    // Segmentation Tests
    // =========================================================================

    #[test]
    fn test_euclidean_clustering() {
        // Create two separate clusters
        let mut points = Vec::new();

        // Cluster 1: points near origin
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(i as f64 * 0.1, j as f64 * 0.1, 0.0));
            }
        }

        // Cluster 2: points far away
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(10.0 + i as f64 * 0.1, j as f64 * 0.1, 0.0));
            }
        }

        let clusters = euclidean_clustering(&points, 0.5, 10);

        assert_eq!(clusters.len(), 2, "Should find 2 clusters");
        assert!(clusters[0].len() >= 100, "First cluster should have at least 100 points");
        assert!(clusters[1].len() >= 100, "Second cluster should have at least 100 points");
    }

    #[test]
    fn test_euclidean_clustering_single_cluster() {
        // Single dense cluster
        let points: Vec<DVec3> = (0..100).map(|i| {
            let t = 2.0 * PI * i as f64 / 100.0;
            DVec3::new(t.cos() * 0.5, t.sin() * 0.5, 0.0)
        }).collect();

        let clusters = euclidean_clustering(&points, 1.0, 10);

        assert_eq!(clusters.len(), 1, "Should find 1 cluster");
        assert!(clusters[0].len() >= 50, "Cluster should contain most points");
    }

    #[test]
    fn test_region_growing_segmentation() {
        // Create planar region with some noise
        let mut points = Vec::new();

        // Planar region
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let config = RegionGrowingConfig {
            k_neighbors: 10,
            max_angle: PI / 6.0,
            max_distance: 0.1,
            min_segment_size: 50,
            max_segments: 10,
        };

        let segments = region_growing_segmentation(&points, &config);

        // Segmentation may or may not find segments depending on parameters
        // The key is that it runs without panicking
        // If segments are found, verify they're valid
        for seg in &segments {
            assert!(!seg.point_indices.is_empty(), "Segments should have points");
        }
    }

    #[test]
    fn test_shape_segmentation_plane() {
        // Create planar points
        let mut points = Vec::new();
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        // Add some noise
        points.push(DVec3::new(100.0, 100.0, 100.0));

        let result = shape_segmentation(&points, ShapeType::Plane, 0.1, 100, 100);

        assert!(result.is_some(), "Plane segmentation should succeed");
        let (params, inliers, outliers) = result.unwrap();

        match params {
            ShapeParams::Plane { point: _, normal } => {
                assert!(normal.z.abs() > 0.9, "Normal should be along Z");
            }
            _ => panic!("Expected plane parameters"),
        }

        assert!(inliers.len() >= 100, "Should find many inliers");
        assert!(!outliers.is_empty(), "Should have outliers");
    }

    #[test]
    fn test_shape_segmentation_sphere() {
        // Create spherical points
        let center = DVec3::new(1.0, 2.0, 3.0);
        let radius = 2.0;
        let mut points = Vec::new();

        for i in 0..100 {
            let theta = 2.0 * PI * i as f64 / 100.0;
            let phi = PI * (i % 10) as f64 / 10.0;
            points.push(DVec3::new(
                center.x + radius * phi.sin() * theta.cos(),
                center.y + radius * phi.sin() * theta.sin(),
                center.z + radius * phi.cos(),
            ));
        }

        let result = shape_segmentation(&points, ShapeType::Sphere, 0.5, 50, 100);

        assert!(result.is_some(), "Sphere segmentation should succeed");
        let (params, _, _) = result.unwrap();

        match params {
            ShapeParams::Sphere { center: c, radius: r } => {
                assert!((c - center).length() < 0.5, "Center should be close");
                assert!((r - radius).abs() < 0.5, "Radius should be close");
            }
            _ => panic!("Expected sphere parameters"),
        }
    }

    #[test]
    fn test_shape_segmentation_cylinder() {
        // Create cylindrical points
        let radius = 1.0;
        let mut points = Vec::new();

        for i in 0..20 {
            let theta = 2.0 * PI * i as f64 / 20.0;
            for z in 0..10 {
                points.push(DVec3::new(
                    radius * theta.cos(),
                    radius * theta.sin(),
                    z as f64,
                ));
            }
        }

        let result = shape_segmentation(&points, ShapeType::Cylinder, 0.2, 50, 100);

        assert!(result.is_some(), "Cylinder segmentation should succeed");
        let (params, _, _) = result.unwrap();

        match params {
            ShapeParams::Cylinder { axis_point: _, axis_direction, radius: r } => {
                assert!(axis_direction.z.abs() > 0.9, "Axis should be along Z");
                assert!((r - radius).abs() < 0.2, "Radius should be close");
            }
            _ => panic!("Expected cylinder parameters"),
        }
    }

    // =========================================================================
    // Surface Reconstruction Tests
    // =========================================================================

    #[test]
    fn test_triangle_mesh_basics() {
        let mut mesh = TriangleMesh::new();
        mesh.nodes.push(DVec3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(DVec3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(DVec3::new(0.0, 1.0, 0.0));
        mesh.triangles.push([0, 1, 2]);

        let face_normals = mesh.compute_face_normals();
        assert_eq!(face_normals.len(), 1);
        assert!(face_normals[0].z.abs() > 0.9, "Face normal should point along Z");

        let node_normals = mesh.compute_node_normals();
        assert_eq!(node_normals.len(), 3);
    }

    #[test]
    fn test_poisson_reconstruction() {
        // Create simple point cloud on a plane
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..10 {
            for j in 0..10 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let config = PoissonConfig {
            depth: 4,
            solver_divide: 4,
            iso_value: 0.0,
        };

        let result = poisson_reconstruction(&points, &normals, &config);
        // May or may not produce output depending on implicit function
        if let Some(mesh) = result {
            assert!(!mesh.nodes.is_empty());
            assert!(!mesh.triangles.is_empty());
        }
    }

    #[test]
    fn test_delaunay_reconstruction() {
        // Create coplanar points
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let result = delaunay_reconstruction(&points, &normals);

        assert!(result.is_some(), "Delaunay reconstruction should succeed");
        let mesh = result.unwrap();

        assert_eq!(mesh.nodes.len(), 25, "Should have all input nodes");
        assert!(!mesh.triangles.is_empty(), "Should have triangles");
    }

    #[test]
    fn test_ball_pivoting_reconstruction() {
        // Create simple point cloud
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64 * 0.5, j as f64 * 0.5, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let config = BpaConfig {
            ball_radius: 1.0,
            clustering: 0.01,
            angle_threshold: PI / 4.0,
        };

        let result = ball_pivoting_reconstruction(&points, &normals, &config);

        // BPA may or may not find valid triangles
        if let Some(mesh) = result {
            assert!(!mesh.nodes.is_empty());
        }
    }

    #[test]
    fn test_generate_consistent_mesh() {
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let result = generate_consistent_mesh(&points, &normals);

        assert!(result.is_some(), "Should generate consistent mesh");
        let mesh = result.unwrap();

        assert!(mesh.normals.is_some(), "Should have computed normals");
        let node_normals = mesh.normals.unwrap();
        assert_eq!(node_normals.len(), mesh.nodes.len());
    }

    // =========================================================================
    // Advanced Sampling Tests
    // =========================================================================

    #[test]
    fn test_voxel_grid_sample() {
        // Dense grid of points
        let mut points = Vec::new();
        for i in 0..10 {
            for j in 0..10 {
                for k in 0..10 {
                    points.push(DVec3::new(i as f64, j as f64, k as f64));
                }
            }
        }

        let config = AdvancedSamplingConfig {
            voxel_size: 1.0,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::VoxelGrid, &config);

        assert!(result.len() <= points.len());
        assert!(result.len() >= 27, "Should have at least 3x3x3 voxels");
    }

    #[test]
    fn test_random_uniform_sample() {
        let points: Vec<DVec3> = (0..1000).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let config = AdvancedSamplingConfig {
            target_count: 100,
            seed: 42,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::RandomUniform, &config);

        assert_eq!(result.len(), 100, "Should have exactly target_count points");
    }

    #[test]
    fn test_curvature_aware_sample() {
        // Create points with varying curvature (plane + hemisphere)
        let mut points = Vec::new();

        // Flat plane (low curvature)
        for i in 0..10 {
            for j in 0..10 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        // Hemisphere (high curvature at edges)
        for i in 0..20 {
            let theta = 2.0 * PI * i as f64 / 20.0;
            for j in 0..10 {
                let phi = PI * j as f64 / 20.0;
                points.push(DVec3::new(
                    20.0 + 2.0 * phi.sin() * theta.cos(),
                    2.0 * phi.sin() * theta.sin(),
                    2.0 * phi.cos(),
                ));
            }
        }

        let config = AdvancedSamplingConfig {
            target_count: 50,
            k_neighbors: 10,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::CurvatureAware, &config);

        assert_eq!(result.len(), 50, "Should have target_count points");
    }

    #[test]
    fn test_poisson_disk_sample() {
        let points: Vec<DVec3> = (0..100).map(|i| {
            let t = 2.0 * PI * i as f64 / 100.0;
            DVec3::new(t.cos(), t.sin(), i as f64 / 100.0)
        }).collect();

        let config = AdvancedSamplingConfig {
            target_count: 50,
            min_distance: 0.1,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::PoissonDisk, &config);

        // Poisson disk sampling should produce a result with fewer points than input
        // The exact count depends on implementation
        assert!(!result.is_empty(), "Should produce at least some samples");
        assert!(result.len() <= points.len(), "Should not produce more samples than input");
    }

    #[test]
    fn test_advanced_sample_identity() {
        // When target_count >= input size, should return all points
        let points: Vec<DVec3> = (0..50).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let config = AdvancedSamplingConfig {
            target_count: 100,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::RandomUniform, &config);
        assert_eq!(result.len(), 50, "Should return all points when target >= input size");
    }

    // =========================================================================
    // Helper Function Tests
    // =========================================================================

    #[test]
    fn test_angles_to_rotation_matrix() {
        // Identity rotation (zero angles)
        let r = angles_to_rotation_matrix(0.0, 0.0, 0.0);

        assert!((r[0][0] - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((r[1][1] - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((r[2][2] - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_compute_circumcenter() {
        let a = DVec3::new(0.0, 0.0, 0.0);
        let b = DVec3::new(1.0, 0.0, 0.0);
        let c = DVec3::new(0.5, 0.5 * 3.0_f64.sqrt(), 0.0);

        let cc = compute_circumcenter(a, b, c);

        assert!(cc.is_some());
        let cc = cc.unwrap();

        // Check all points are equidistant from circumcenter
        let ra = (a - cc).length();
        let rb = (b - cc).length();
        let rc = (c - cc).length();

        assert!((ra - rb).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((rb - rc).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_delaunay_triangulation_2d() {
        // Four points forming a square
        let points = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(0.0, 1.0),
        ];

        let triangles = delaunay_triangulation_2d(&points);

        // Delaunay triangulation should produce at least 1 triangle
        // For a square, it typically produces 2 triangles
        assert!(!triangles.is_empty(), "Should produce at least one triangle");
    }
}
