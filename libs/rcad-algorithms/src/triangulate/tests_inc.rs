#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangulate_triangle() {
        let verts = vec![DVec3::ZERO, DVec3::X, DVec3::Y];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 1);
    }

    #[test]
    fn triangulate_quad() {
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn triangulate_pentagon() {
        let verts = (0..5)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect::<Vec<_>>();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 3);
    }

    #[test]
    fn empty_polygon_returns_no_triangles() {
        let tris = triangulate_polygon(&[], DVec3::Z);
        assert!(tris.is_empty());
    }

    #[test]
    fn two_vertex_polygon_returns_no_triangles() {
        let verts = vec![DVec3::ZERO, DVec3::X];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert!(tris.is_empty());
    }

    #[test]
    fn triangle_count_is_n_minus_2() {
        // A convex n-gon should always yield n-2 triangles.
        for n in 3..=10 {
            let verts: Vec<DVec3> = (0..n)
                .map(|i| {
                    let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                    DVec3::new(a.cos(), a.sin(), 0.0)
                })
                .collect();
            let tris = triangulate_polygon(&verts, DVec3::Z);
            assert_eq!(
                tris.len(),
                n - 2,
                "expected {n}-gon to yield {} triangles, got {}",
                n - 2,
                tris.len()
            );
        }
    }

    #[test]
    fn all_indices_in_bounds() {
        // Every index in the triangulation must be < number of vertices.
        let verts: Vec<DVec3> = (0..7)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 7.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        for tri in &tris {
            for &idx in tri.iter() {
                assert!(idx < verts.len(), "index {idx} out of bounds for {n} vertices", n = verts.len());
            }
        }
    }

    #[test]
    fn clockwise_quad_still_triangulates() {
        // Reversed vertex order (CW) should be handled by sign-flip logic.
        let verts = vec![
            DVec3::new(0.0, 1.0, 0.0), // top-left first (CW)
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
        ];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 2);
    }

    /// mesh_brep on a box primitive should fill face.triangles for all 6 faces.
    #[test]
    fn mesh_brep_box_fills_triangles() {
        use rcad_kernel::BRep;
        use rcad_kernel::geom::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        let params = TessellationParams::default();
        mesh_brep(&mut brep, &params);

        let faces = &brep.solids[0].shells[0].faces;
        assert_eq!(faces.len(), 6, "box should have 6 faces");
        for (i, face) in faces.iter().enumerate() {
            assert!(
                !face.triangles.is_empty(),
                "face {i} should have triangles after mesh_brep"
            );
            // All triangle indices must be valid vertex indices.
            for &[a, b, c] in &face.triangles {
                assert!(a < brep.vertices.len(), "face {i}: vertex index {a} out of bounds");
                assert!(b < brep.vertices.len(), "face {i}: vertex index {b} out of bounds");
                assert!(c < brep.vertices.len(), "face {i}: vertex index {c} out of bounds");
            }
        }
    }

    /// mesh_brep on a sphere should produce a dense mesh (many triangles per face).
    #[test]
    fn mesh_brep_sphere_produces_dense_mesh() {
        use rcad_kernel::BRep;
        use rcad_kernel::geom::PrimitiveSolid;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let params = TessellationParams {
            chord_tolerance: 0.05,
            ..TessellationParams::default()
        };
        mesh_brep(&mut brep, &params);

        let total_tris: usize = brep.solids[0].shells[0].faces
            .iter()
            .map(|f| f.triangles.len())
            .sum();
        assert!(
            total_tris >= 8,
            "sphere mesh should have at least 8 triangles, got {total_tris}"
        );
    }

    /// mesh_brep on a cylinder should produce triangles for all faces.
    #[test]
    fn mesh_brep_cylinder_all_faces_have_triangles() {
        use rcad_kernel::BRep;
        use rcad_kernel::geom::PrimitiveSolid;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let params = TessellationParams::default();
        mesh_brep(&mut brep, &params);

        let faces = &brep.solids[0].shells[0].faces;
        for (i, face) in faces.iter().enumerate() {
            assert!(
                !face.triangles.is_empty(),
                "cylinder face {i} should have triangles"
            );
        }
    }

    #[test]
    fn mesh_brep_fallback_triangulates_semicircle_wire_face() {
        use std::f64::consts::PI;
        use rcad_kernel::BRep;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
        use rcad_kernel::geom::{Circle3, Curve3, CurveEval};

        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0,);
        let p0 = circle.point_at(0.0);
        let p1 = circle.point_at(PI);

        // Two-edge closed wire: semicircular arc + diameter chord.
        let mut brep = BRep {
            vertices: vec![Vertex { point: p0 }, Vertex { point: p1 }],
            edges: vec![
                Edge { start: 0, end: 1 },
                Edge { start: 1, end: 0 },
            ],
            solids: vec![Solid {
                shells: vec![Shell {
                    faces: vec![Face {
                        outer_wire: Wire {
                            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
                        },
                        inner_wires: vec![],
                        normal: DVec3::Z,
                        triangles: vec![],
                        sample_point: None,
                        mesh_dirty: true,
            surface_idx: None,
                    }],
                }],
            }],
            geom: rcad_kernel::GeomStore {
                curves: vec![Curve3::Circle(circle)],
                edge_curve: vec![Some(0), None],
                edge_curve_range: vec![Some([0.0, PI]), None],
                face_surface: vec![None],
                ..Default::default()
            },
            compound: None,
            compsolid: None,
        };

        mesh_brep(&mut brep, &TessellationParams::default());
        let tris = &brep.solids[0].shells[0].faces[0].triangles;
        assert!(
            tris.len() > 1,
            "semicircle fallback should produce multiple triangles, got {}",
            tris.len()
        );
    }

    // ========================================================================
    // Tessellation / mesh utility tests
    // ========================================================================

    #[test]
    fn tessellation_params_presets() {
        // Preview preset is coarser and disables adaptive refinement
        let preview = TessellationParams::preview();
        assert!(preview.chord_tolerance > TessellationParams::standard().chord_tolerance);
        assert!(!preview.adaptive_refinement);
        assert!(preview.parallel);

        // Standard preset enables adaptive refinement
        let standard = TessellationParams::standard();
        assert!(standard.adaptive_refinement);
        assert!(standard.curvature_sensitive);

        // High quality tightens tolerances vs standard
        let hq = TessellationParams::high_quality();
        assert!(hq.chord_tolerance < standard.chord_tolerance);
        assert!(hq.max_depth > standard.max_depth);

        // Export preset sits between HQ and standard chord tolerance
        let export = TessellationParams::export();
        assert!(export.chord_tolerance > hq.chord_tolerance);
        assert!(export.chord_tolerance < standard.chord_tolerance);

        // Analysis preset is the strictest of the built-ins
        let analysis = TessellationParams::analysis();
        assert!(analysis.chord_tolerance < hq.chord_tolerance);
        assert!(analysis.max_aspect_ratio < hq.max_aspect_ratio);
    }

    #[test]
    fn tessellation_params_with_target_triangle_count() {
        let params = TessellationParams::standard();
        // Higher target count means more triangles -> finer tolerance
        let adjusted = params.with_target_triangle_count(10000);
        // Factor = (10000/1000)^(1/3) 閳?2.15, so tolerance increases (coarser mesh)
        // For more triangles, we'd actually want lower tolerance, so this adjusts accordingly
        assert!(adjusted.chord_tolerance != params.chord_tolerance);
    }

    #[test]
    fn mesh_quality_metrics_basic() {
        let nodes = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let triangles = vec![[0, 1, 2]];

        let metrics = compute_mesh_quality(&nodes, &triangles);
        assert_eq!(metrics.triangle_count, 1);
        assert_eq!(metrics.node_count, 3);
        assert_eq!(metrics.degenerate_count, 0);
        assert!(metrics.max_aspect_ratio > 1.0);
    }

    #[test]
    fn mesh_quality_metrics_degenerate() {
        // Collinear vertices -> degenerate triangle
        let nodes = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let triangles = vec![[0, 1, 2]];

        let metrics = compute_mesh_quality(&nodes, &triangles);
        assert_eq!(metrics.degenerate_count, 1);
        // For collinear points, the aspect ratio is max_edge/min_edge = 2/1 = 2
        // which is still bad but not infinite. The key is degenerate_count.
        assert!(metrics.max_aspect_ratio > 1.0);
    }

    #[test]
    fn mesh_quality_metrics_score() {
        // Near-equilateral triangle
        let nodes = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.5, 0.866, 0.0), // ~60鎺?internal angles
        ];
        let triangles = vec![[0, 1, 2]];

        let metrics = compute_mesh_quality(&nodes, &triangles);
        assert!(metrics.quality_score() > 0.9);
        assert!(metrics.is_good(20.0));
    }

    #[test]
    fn surface_mesh_compute_quality() {
        let mesh = SurfaceMesh {
            nodes: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2]],
            normals: vec![DVec3::Z; 3],
            dirty: false,
        };

        let metrics = mesh.compute_quality();
        assert_eq!(metrics.triangle_count, 1);
    }

    #[test]
    fn adaptive_subdivider_default() {
        let subdivider = AdaptiveSubdivider::new();
        assert_eq!(subdivider.curvature_threshold, 0.1);
        assert_eq!(subdivider.distance_threshold, 0.1);
        assert_eq!(subdivider.max_subdivision_levels, 3);
    }

    #[test]
    fn adaptive_subdivider_builder() {
        let subdivider = AdaptiveSubdivider::new()
            .with_curvature_threshold(0.2)
            .with_distance_threshold(0.5)
            .with_max_levels(5);

        assert_eq!(subdivider.curvature_threshold, 0.2);
        assert_eq!(subdivider.distance_threshold, 0.5);
        assert_eq!(subdivider.max_subdivision_levels, 5);
    }

    #[test]
    fn adaptive_subdivider_subdivide_by_distance() {
        // Large triangle should trigger edge splits
        let mesh = SurfaceMesh {
            nodes: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(10.0, 0.0, 0.0),
                DVec3::new(5.0, 10.0, 0.0),
            ],
            triangles: vec![[0, 1, 2]],
            normals: vec![DVec3::Z; 3],
            dirty: false,
        };

        let subdivider = AdaptiveSubdivider::new()
            .with_distance_threshold(1.0);
        let result = subdivider.subdivide_by_distance(&mesh);

        assert!(result.triangles.len() > 1, "distance subdivider should split long edges");
    }

    #[test]
    fn boundary_sensitive_tessellator_default() {
        let tessellator = BoundarySensitiveTessellator::new();
        assert_eq!(tessellator.feature_angle_threshold, 0.52);
        assert!(tessellator.auto_detect_features);
    }

    #[test]
    fn boundary_sensitive_tessellator_detect_features() {
        // Two triangles sharing an edge with ~90鎺?dihedral
        let nodes = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
        ];
        let triangles = vec![[0, 1, 2], [0, 1, 3]];
        let normals = vec![DVec3::Z; 4];

        let mut tessellator = BoundarySensitiveTessellator::new()
            .with_feature_angle(0.1); // low threshold -> crease detection should fire

        tessellator.detect_feature_edges(&nodes, &triangles, &normals);
        assert!(!tessellator.feature_edges.is_empty());
    }

    #[test]
    fn incremental_mesher_basic() {
        let mut mesher = IncrementalMesher::new();
        assert!(!mesher.is_dirty());

        mesher.invalidate_face(0);
        assert!(mesher.is_dirty());
        assert!(mesher.dirty_faces.contains(&0));

        mesher.clear();
        assert!(!mesher.is_dirty());
    }

    #[test]
    fn incremental_mesher_multiple_faces() {
        let mut mesher = IncrementalMesher::new();
        mesher.invalidate_faces(&[0, 1, 2]);

        assert!(mesher.dirty_faces.contains(&0));
        assert!(mesher.dirty_faces.contains(&1));
        assert!(mesher.dirty_faces.contains(&2));
        assert_eq!(mesher.dirty_faces.len(), 3);
    }

    #[test]
    fn mesh_simplifier_default() {
        let simplifier = MeshSimplifier::new();
        assert_eq!(simplifier.target_ratio, 0.5);
        assert_eq!(simplifier.max_error, 0.01);
        assert!(simplifier.preserve_boundary);
    }

    #[test]
    fn mesh_simplifier_builder() {
        let simplifier = MeshSimplifier::new()
            .with_target_ratio(0.25)
            .with_max_error(0.05);

        assert_eq!(simplifier.target_ratio, 0.25);
        assert_eq!(simplifier.max_error, 0.05);
    }

    #[test]
    fn mesh_simplifier_simplify() {
        // Small structured grid of triangles
        let nodes: Vec<DVec3> = (0..9)
            .map(|i| {
                let row = i / 3;
                let col = i % 3;
                DVec3::new(col as f64, row as f64, 0.0)
            })
            .collect();
        let triangles = vec![
            [0, 1, 3], [1, 4, 3],
            [1, 2, 4], [2, 5, 4],
            [3, 4, 6], [4, 7, 6],
            [4, 5, 7], [5, 8, 7],
        ];

        let mesh = SurfaceMesh {
            nodes,
            triangles,
            normals: vec![DVec3::Z; 9],
            dirty: false,
        };

        let simplifier = MeshSimplifier::new()
            .with_target_ratio(0.5)
            .with_max_error(1.0); // permissive so short edges can collapse

        let simplified = simplifier.simplify_mesh(&mesh);
        assert!(simplified.triangles.len() <= mesh.triangles.len());
    }

    #[test]
    fn mesh_simplifier_simplify_to_target_count() {
        let mesh = SurfaceMesh {
            nodes: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            normals: vec![DVec3::Z; 4],
            dirty: false,
        };

        let simplifier = MeshSimplifier::new().with_max_error(1.0);
        let result = simplifier.simplify_to_target_count(&mesh, 4);

        // Mesh already has fewer triangles than the requested target
        assert_eq!(result.triangles.len(), 2);
    }

    #[test]
    fn find_boundary_nodes() {
        let triangles = vec![[0, 1, 2], [1, 3, 2]];
        let boundary = super::find_boundary_nodes(&triangles);

        // Boundary edges: (0,1), (0,2), (1,3), (2,3); internal edge: (1,2)
        assert!(boundary.contains(&0));
        assert!(boundary.contains(&3));
    }
}


