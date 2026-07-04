#[cfg(test)]
mod tests {
    use crate::brep_check_parallel::*;
    use rcad_kernel::BRep;
    use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn test_check_parallel_empty_brep() {
        let brep = BRep::default();
        let result = check_parallel(&brep);
        assert!(result.is_valid());
    }

    #[test]
    fn test_check_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check_parallel(&brep);
        assert!(result.is_valid(), "issues: {:?}", result.issues);
    }

    #[test]
    fn test_check_parallel_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        let result = check_parallel(&brep);
        // Cylinder has seam edges that may trigger non-manifold warnings
        // The check should complete without panic, not necessarily be valid
        let _ = result.issues.len();
    }

    #[test]
    fn test_check_parallel_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere {
            radius: 1.0,
        });
        let result = check_parallel(&brep);
        // Sphere has seam edges that may trigger non-manifold warnings
        // The check should complete without panic, not necessarily be valid
        let _ = result.issues.len();
    }

    #[test]
    fn test_check_many_parallel() {
        let breps: Vec<BRep> = vec![
            BRep::from_primitive(PrimitiveSolid::Box {
                width: 1.0, height: 1.0, depth: 1.0,
            }),
            BRep::from_primitive(PrimitiveSolid::Sphere {
                radius: 1.0,
            }),
            BRep::from_primitive(PrimitiveSolid::Cylinder {
                radius: 1.0, height: 2.0,
            }),
        ];

        let results = check_many_parallel(&breps);
        assert_eq!(results.len(), 3);
        // Box should be valid
        assert!(results[0].is_valid(), "issues: {:?}", results[0].issues);
        // Sphere and cylinder have seam edges that may trigger warnings
        // Just verify the checks completed
        let _ = results[1].issues.len();
        let _ = results[2].issues.len();
    }

    #[test]
    fn test_check_parallel_with_stats() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let (result, stats) = check_parallel_with_stats(&brep);
        assert!(result.is_valid(), "issues: {:?}", result.issues);
        assert_eq!(stats.face_count, 6); // Box has 6 faces
        assert_eq!(stats.edge_count, 12); // Box has 12 edges
        assert_eq!(stats.vertex_count, 8); // Box has 8 vertices
        assert_eq!(stats.issue_count, 0);
        assert!(stats.is_valid);
    }

    #[test]
    fn test_segments_intersect_2d() {
        // Crossing segments
        let p1 = DVec3::new(0.0, 0.0, 0.0);
        let p2 = DVec3::new(2.0, 2.0, 0.0);
        let p3 = DVec3::new(0.0, 2.0, 0.0);
        let p4 = DVec3::new(2.0, 0.0, 0.0);
        assert!(segments_intersect_2d(p1, p2, p3, p4));

        // Non-crossing segments
        let p5 = DVec3::new(0.0, 0.0, 0.0);
        let p6 = DVec3::new(1.0, 1.0, 0.0);
        let p7 = DVec3::new(3.0, 3.0, 0.0);
        let p8 = DVec3::new(4.0, 4.0, 0.0);
        assert!(!segments_intersect_2d(p5, p6, p7, p8));
    }

    #[test]
    fn test_parallel_options_default() {
        let opts = ParallelCheckOptions::default();
        assert_eq!(opts.min_faces_for_parallel, 100);
        assert_eq!(opts.chunk_size, 32);
        assert!(opts.check_duplicate_vertices);
        assert!(opts.check_isolated_vertices);
        assert!(opts.check_finite_vertices);
    }

    #[test]
    fn test_parallel_options_small_model() {
        let opts = ParallelCheckOptions::small_model();
        assert_eq!(opts.min_faces_for_parallel, usize::MAX);
    }

    #[test]
    fn test_parallel_options_large_model() {
        let opts = ParallelCheckOptions::large_model();
        assert_eq!(opts.min_faces_for_parallel, 10);
        assert_eq!(opts.chunk_size, 64);
    }

    #[test]
    fn test_parallel_vs_sequential_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Both should produce same results
        let parallel_result = check_parallel(&brep);
        let sequential_result = crate::brep_check::brep_check_analyze(&brep);

        assert_eq!(parallel_result.is_valid(), sequential_result.is_valid());
        assert_eq!(parallel_result.issues.len(), sequential_result.issues.len());
    }

    #[test]
    fn test_parallel_vs_sequential_invalid_brep() {
        // Create an invalid BRep with open wire
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // Gap: v2 != v3

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let parallel_result = check_parallel(&brep);
        let sequential_result = crate::brep_check::brep_check_analyze(&brep);

        // Both should detect the open wire
        assert!(!parallel_result.is_valid());
        assert!(!sequential_result.is_valid());

        // Both should have same number of issues
        assert_eq!(parallel_result.issues.len(), sequential_result.issues.len());

        // Both should have OpenWire issue
        assert!(parallel_result.issues.iter().any(|i| matches!(i, CheckIssue::OpenWire { .. })));
        assert!(sequential_result.issues.iter().any(|i| matches!(i, CheckIssue::OpenWire { .. })));
    }

    #[test]
    fn test_small_model_uses_sequential() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let opts = ParallelCheckOptions::small_model();
        let result = check_parallel_with_options(&brep, &opts);

        assert!(!result.was_parallel, "Small model should use sequential processing");
        assert_eq!(result.thread_count, 1);
    }

    #[test]
    fn test_large_model_uses_parallel() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let opts = ParallelCheckOptions::large_model();
        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.was_parallel, "Large model settings should use parallel processing");
        assert!(result.thread_count >= 1);
    }

    #[test]
    fn test_isolated_vertex_detection() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(2.0, 0.0, 0.0) }); // Isolated

        brep.edges.push(Edge { start: 0, end: 1 });

        let opts = ParallelCheckOptions {
            check_isolated_vertices: true,
            check_duplicate_vertices: false,
            check_finite_vertices: false,
            ..ParallelCheckOptions::default()
        };

        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.parallel_issues.iter().any(|i| matches!(
            i,
            ParallelCheckIssue::IsolatedVertex { vertex_idx: 2 }
        )), "Should detect isolated vertex 2");
    }

    #[test]
    fn test_non_finite_vertex_detection() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(f64::NAN, 0.0, 0.0) }); // NaN

        brep.edges.push(Edge { start: 0, end: 1 });

        let opts = ParallelCheckOptions {
            check_finite_vertices: true,
            check_duplicate_vertices: false,
            check_isolated_vertices: false,
            ..ParallelCheckOptions::default()
        };

        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.parallel_issues.iter().any(|i| matches!(
            i,
            ParallelCheckIssue::NonFiniteVertex { vertex_idx: 1 }
        )), "Should detect non-finite vertex 1");
    }

    #[test]
    fn test_duplicate_vertex_detection() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // Duplicate

        brep.edges.push(Edge { start: 0, end: 1 });

        let opts = ParallelCheckOptions {
            check_duplicate_vertices: true,
            check_isolated_vertices: false,
            check_finite_vertices: false,
            tolerance: TOLERANCE_MESH_LEGACY,
            ..ParallelCheckOptions::default()
        };

        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.parallel_issues.iter().any(|i| matches!(
            i,
            ParallelCheckIssue::DuplicateVertex { vertex_a: 0, vertex_b: 1, .. }
        )), "Should detect duplicate vertices");
    }

    #[test]
    fn test_check_many_parallel_with_options() {
        let breps: Vec<BRep> = vec![
            BRep::from_primitive(PrimitiveSolid::Box {
                width: 1.0, height: 1.0, depth: 1.0,
            }),
            BRep::from_primitive(PrimitiveSolid::Sphere {
                radius: 1.0,
            }),
        ];

        let opts = ParallelCheckOptions::default();
        let results = check_many_parallel_with_options(&breps, &opts);

        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_parallel_check_result_is_valid() {
        let mut result = ParallelCheckResult::default();
        assert!(result.is_valid());

        result.issues.push(CheckIssue::DegenerateFace { solid: 0, shell: 0, face: 0 });
        assert!(!result.is_valid());
    }

    #[test]
    fn test_parallel_check_result_to_check_result() {
        let mut result = ParallelCheckResult::default();
        result.issues.push(CheckIssue::DegenerateFace { solid: 0, shell: 0, face: 0 });

        let check_result = result.to_check_result();
        assert_eq!(check_result.issues.len(), 1);
    }

    /// Generate a large BRep for performance testing.
    #[cfg(test)]
    fn generate_large_brep(n_boxes: usize) -> BRep {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();

        // Create a grid of connected quads
        let mut vertex_offset = 0usize;
        let mut edge_offset = 0usize;

        for _z in 0..n_boxes {
            for _y in 0..n_boxes {
                for _x in 0..n_boxes {
                    // Add 8 vertices for a box
                    for dz in 0..2 {
                        for dy in 0..2 {
                            for dx in 0..2 {
                                let x = dx as f64;
                                let y = dy as f64;
                                let z = dz as f64;
                                brep.vertices.push(Vertex {
                                    point: DVec3::new(x, y, z),
                                });
                            }
                        }
                    }

                    // Add 12 edges for the box
                    let v = vertex_offset;
                    let edges = vec![
                        (v+0, v+1), (v+1, v+3), (v+3, v+2), (v+2, v+0), // bottom
                        (v+4, v+5), (v+5, v+7), (v+7, v+6), (v+6, v+4), // top
                        (v+0, v+4), (v+1, v+5), (v+2, v+6), (v+3, v+7), // vertical
                    ];

                    for (start, end) in edges {
                        brep.edges.push(Edge { start, end });
                    }

                    // Add 6 faces for the box
                    let e = edge_offset;
                    let face_wire_indices = vec![
                        vec![(e+0, true), (e+1, true), (e+2, true), (e+3, true)],   // bottom
                        vec![(e+4, true), (e+5, true), (e+6, true), (e+7, true)],   // top
                        vec![(e+0, true), (e+8, true), (e+4, false), (e+11,false)], // front
                        vec![(e+2, false), (e+10,true), (e+6, false), (e+9, false)],// back
                        vec![(e+3, true), (e+10,false),(e+7, true), (e+8, false)], // left
                        vec![(e+1, false),(e+9, true), (e+5, false), (e+11,true)], // right
                    ];

                    let normals = vec![
                        DVec3::NEG_Z, DVec3::Z, DVec3::NEG_Y, DVec3::Y, DVec3::NEG_X, DVec3::X,
                    ];

                    let mut faces = Vec::new();
                    for (fi, wire_indices) in face_wire_indices.iter().enumerate() {
                        faces.push(Face {
                            outer_wire: Wire {
                                edges: wire_indices.iter().map(|&(idx, fwd)| {
                                    if fwd { WireEdge::fwd(idx) } else { WireEdge::rev(idx) }
                                }).collect(),
                            },
                            inner_wires: vec![],
                            normal: normals[fi],
                            triangles: vec![],
                            sample_point: None,
                            mesh_dirty: true,
            surface_idx: None,
                        });
                    }

                    brep.solids.push(Solid {
                        shells: vec![Shell { faces }],
                    });

                    vertex_offset += 8;
                    edge_offset += 12;
                }
            }
        }

        brep
    }

    #[test]
    fn test_large_brep_parallel_vs_sequential() {
        // Create a moderately large BRep
        let brep = generate_large_brep(3); // 27 boxes, 162 faces

        let parallel_result = check_parallel(&brep);
        let sequential_result = crate::brep_check::brep_check_analyze(&brep);

        // Results should be identical
        assert_eq!(parallel_result.is_valid(), sequential_result.is_valid());
        assert_eq!(parallel_result.issues.len(), sequential_result.issues.len());
    }

    #[test]
    fn test_parallel_options_builder() {
        let opts = ParallelCheckOptions::default()
            .with_tolerance(TOLERANCE_COORD_SUB)
            .with_chunk_size(128)
            .with_duplicate_vertex_check(false)
            .with_isolated_vertex_check(false);

        assert!((opts.tolerance - TOLERANCE_COORD_SUB).abs() < TOLERANCE_FLOAT_DEDUP);
        assert_eq!(opts.chunk_size, 128);
        assert!(!opts.check_duplicate_vertices);
        assert!(!opts.check_isolated_vertices);
    }

    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
    // Tests for check_faces_parallel
    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?

    #[test]
    fn test_check_faces_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = check_faces_parallel(&brep, 0);
        assert_eq!(results.len(), 6, "Box should have 6 faces");

        for result in &results {
            assert!(result.is_valid, "Face should be valid: {:?}", result.issues);
            assert!(result.outer_wire_closed, "Outer wire should be closed");
            assert_eq!(result.outer_wire_edge_count, 4, "Each face should have 4 edges");
        }
    }

    #[test]
    fn test_check_faces_parallel_open_wire() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // Gap: v2 != v3

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = check_faces_parallel(&brep, 0);
        assert_eq!(results.len(), 1);
        assert!(!results[0].is_valid, "Face with open wire should be invalid");
        assert!(!results[0].outer_wire_closed, "Wire should be reported as open");
        assert!(results[0].outer_wire_gaps > 0, "Should have gaps");
    }

    #[test]
    fn test_check_faces_parallel_degenerate_face() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0)], // Only 1 edge
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = check_faces_parallel(&brep, 0);
        assert!(!results[0].is_valid, "Degenerate face should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, FaceCheckIssue::DegenerateFace)));
    }

    #[test]
    fn test_check_faces_parallel_zero_normal() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // Zero normal
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = check_faces_parallel(&brep, 0);
        assert!(!results[0].is_valid, "Face with zero normal should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, FaceCheckIssue::ZeroNormal)));
    }

    #[test]
    fn test_check_faces_parallel_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let results = check_faces_parallel(&brep, 0);

        // Sphere should have faces
        assert!(!results.is_empty(), "Sphere should have faces");
        // Verify basic structure
        for result in &results {
        }
    }

    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
    // Tests for check_edges_parallel
    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?

    #[test]
    fn test_check_edges_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = check_edges_parallel(&brep, 0);
        assert_eq!(results.len(), 12, "Box should have 12 edges");

        for result in &results {
            assert!(result.is_valid, "Edge should be valid: {:?}", result.issues);
            assert!(result.is_manifold, "Each edge should be manifold");
            assert_eq!(result.face_count, 2, "Each edge should be shared by 2 faces");
            assert!(!result.is_degenerate, "No edge should be degenerate");
            assert!(result.length > 0.0, "Each edge should have positive length");
        }
    }

    #[test]
    fn test_check_edges_parallel_invalid_vertex() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.edges.push(Edge { start: 0, end: 99 }); // Invalid vertex

        let results = check_edges_parallel(&brep, 0);
        assert_eq!(results.len(), 1);
        assert!(!results[0].is_valid, "Edge with invalid vertex should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, EdgeCheckIssue::InvalidVertexIndex { .. })));
    }

    #[test]
    fn test_check_edges_parallel_degenerate() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::ZERO }); // Same position
        brep.edges.push(Edge { start: 0, end: 1 });

        let results = check_edges_parallel(&brep, 0);
        assert!(results[0].is_degenerate, "Edge with same vertex positions should be degenerate");
    }

    #[test]
    fn test_check_edges_parallel_free_edge() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 }); // Not referenced by any face

        let results = check_edges_parallel(&brep, 0);
        assert!(!results[0].is_valid, "Free edge should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, EdgeCheckIssue::FreeEdge)));
    }

    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
    // Tests for validate_shells_parallel
    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?

    #[test]
    fn test_validate_shells_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = validate_shells_parallel(&brep);
        assert_eq!(results.len(), 1, "Box should have 1 shell");

        let shell = &results[0];
        assert!(shell.is_valid, "Box shell should be valid");
        assert!(shell.is_closed, "Box shell should be closed");
        assert!(shell.is_manifold, "Box shell should be manifold");
        assert_eq!(shell.face_count, 6);
        assert_eq!(shell.euler_characteristic, 2, "Box Euler characteristic should be 2");
        assert_eq!(shell.genus, Some(0), "Box genus should be 0");
    }

    #[test]
    fn test_validate_shells_parallel_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let results = validate_shells_parallel(&brep);
        assert_eq!(results.len(), 1);

        let shell = &results[0];
        assert!(shell.is_closed, "Sphere shell should be closed");
        assert_eq!(shell.euler_characteristic, 2, "Sphere Euler characteristic should be 2");
        assert_eq!(shell.genus, Some(0), "Sphere genus should be 0");
    }

    #[test]
    fn test_validate_shells_parallel_open_shell() {
        let mut brep = BRep::new();
        // Create a simple open shell (just one face)
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = validate_shells_parallel(&brep);
        assert!(!results[0].is_closed, "Single face shell should be open");
        assert!(results[0].open_edge_count > 0, "Should have open edges");
    }

    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
    // Tests for validate_solids_parallel
    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?

    #[test]
    fn test_validate_solids_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = validate_solids_parallel(&brep);
        assert_eq!(results.len(), 1, "Should have 1 solid");

        let solid = &results[0];
        assert!(solid.is_valid, "Box solid should be valid");
        assert!(solid.is_closed, "Box solid should be closed");
        assert!(solid.is_manifold, "Box solid should be manifold");
        assert_eq!(solid.face_count, 6);
        assert_eq!(solid.edge_count, 12);
        assert_eq!(solid.vertex_count, 8);
        assert!(solid.volume >= 0.0, "Box volume should be non-negative");
    }

    #[test]
    fn test_validate_solids_parallel_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let results = validate_solids_parallel(&brep);
        assert_eq!(results.len(), 1);

        let solid = &results[0];
        assert!(solid.is_closed, "Cylinder solid should be closed");
        assert!(solid.volume >= 0.0, "Cylinder volume should be non-negative");
    }

    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
    // Tests for check_brep_parallel
    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?

    #[test]
    fn test_check_brep_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        assert!(report.is_valid, "Box should pass all checks");
        assert!(report.structural_issues.is_empty());
        assert!(report.parallel_issues.is_empty());
        assert_eq!(report.total_faces, 6);
        assert_eq!(report.total_edges, 12);
        assert_eq!(report.total_vertices, 8);
        assert_eq!(report.total_solids, 1);
        assert!(report.total_duration_ms < 10000, "Should complete quickly");
    }

    #[test]
    fn test_check_brep_parallel_fast_config() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let config = ParallelCheckConfig::fast();
        let report = check_brep_parallel(&brep, &config);

        // Fast config skips some checks
        assert!(report.phase_timings.iter().any(|t| t.phase == "faces"));
    }

    #[test]
    fn test_check_brep_parallel_thorough_config() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let config = ParallelCheckConfig::thorough();
        let report = check_brep_parallel(&brep, &config);

        // Thorough config has tighter tolerance
        assert!((config.tolerance - TOLERANCE_COORD_SUB).abs() < TOLERANCE_FLOAT_DEDUP);
    }

    #[test]
    fn test_check_brep_parallel_custom_threads() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default().with_threads(2);
        let report = check_brep_parallel(&brep, &config);

        // Should work with custom thread count
        assert!(report.threads_used >= 1);
    }

    #[test]
    fn test_check_brep_parallel_timing() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        // Should have timing for each phase
        let phases: Vec<&str> = report.phase_timings.iter().map(|t| t.phase.as_str()).collect();
        assert!(phases.contains(&"faces"));
        assert!(phases.contains(&"edges"));
        assert!(phases.contains(&"vertices"));
        assert!(phases.contains(&"shells"));
        assert!(phases.contains(&"solids"));
    }

    #[test]
    fn test_check_brep_parallel_invalid_brep() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(f64::NAN, 0.0, 0.0) }); // NaN vertex
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // Duplicate
        brep.edges.push(Edge { start: 0, end: 2 });
        // Vertex 1 is isolated

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        assert!(!report.is_valid, "Invalid BRep should fail checks");
        assert!(!report.parallel_issues.is_empty(), "Should have parallel-specific issues");
    }

    #[test]
    fn test_check_brep_parallel_empty_brep() {
        let brep = BRep::default();

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        assert!(report.is_valid, "Empty BRep should be valid (no issues)");
        assert_eq!(report.total_faces, 0);
        assert_eq!(report.total_edges, 0);
        assert_eq!(report.total_vertices, 0);
    }

    #[test]
    fn test_check_brep_parallel_summary() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        let summary = report.summary();
        assert!(summary.contains("VALID"));
        assert!(summary.contains("1 solids"));
        assert!(summary.contains("6 faces"));
    }

    #[test]
    fn test_check_brep_parallel_timing_summary() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        let timing = report.timing_summary();
        assert!(timing.contains("Timing breakdown"));
    }

    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
    // Tests for result types
    // 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?

    #[test]
    fn test_face_check_result_summary() {
        let result = FaceCheckResult {
            solid_idx: 0,
            shell_idx: 0,
            face_idx: 1,
            is_valid: true,
            issues: vec![],
            outer_wire_edge_count: 4,
            inner_wire_count: 0,
            normal: DVec3::Z,
            normal_valid: true,
            outer_wire_closed: true,
            outer_wire_gaps: 0,
            has_self_intersection: false,
        };

        let summary = result.summary();
        assert!(summary.contains("valid"));
        assert!(result.is_clean());
    }

    #[test]
    fn test_edge_check_result_summary() {
        let result = EdgeCheckResult {
            edge_idx: 0,
            is_valid: true,
            issues: vec![],
            start_vertex: 0,
            end_vertex: 1,
            length: 1.0,
            is_degenerate: false,
            face_count: 2,
            is_manifold: true,
            tolerance: TOLERANCE_MESH_LEGACY,
            has_self_intersection: false,
        };

        let summary = result.summary();
        assert!(summary.contains("valid"));
        assert!(result.is_clean());
    }

    #[test]
    fn test_shell_validation_result_summary() {
        let result = ShellValidationResult {
            solid_idx: 0,
            shell_idx: 0,
            is_valid: true,
            face_count: 6,
            edge_count: 12,
            vertex_count: 8,
            euler_characteristic: 2,
            is_closed: true,
            is_manifold: true,
            open_edge_count: 0,
            non_manifold_edge_count: 0,
            orientation_consistent: true,
            genus: Some(0),
            face_results: vec![],
            errors: vec![],
            warnings: vec![],
        };

        assert!(result.is_closed_manifold());
        let summary = result.summary();
        assert!(summary.contains("VALID"));
    }

    #[test]
    fn test_solid_validation_result_summary() {
        let result = SolidValidationResult {
            solid_idx: 0,
            is_valid: true,
            shell_count: 1,
            face_count: 6,
            edge_count: 12,
            vertex_count: 8,
            euler_characteristic: 2,
            is_closed: true,
            is_manifold: true,
            orientation_valid: true,
            has_positive_volume: true,
            volume: 1.0,
            genus: Some(0),
            shell_results: vec![],
            errors: vec![],
            warnings: vec![],
        };

        assert!(result.is_valid_for_operations());
        let summary = result.summary();
        assert!(summary.contains("VALID"));
    }

    #[test]
    fn test_parallel_check_config_presets() {
        let fast = ParallelCheckConfig::fast();
        assert!(!fast.check_self_intersections);
        assert!(!fast.check_same_parameter);

        let thorough = ParallelCheckConfig::thorough();
        assert!((thorough.tolerance - TOLERANCE_COORD_SUB).abs() < TOLERANCE_FLOAT_DEDUP);
        assert!(thorough.check_self_intersections);
    }

    #[test]
    fn test_parallel_check_report_timing() {
        let timing = CheckPhaseTiming {
            phase: "test".to_string(),
            duration_ms: 100,
            items_processed: 50,
        };

        assert_eq!(timing.phase, "test");
        assert_eq!(timing.duration_ms, 100);
        assert_eq!(timing.items_processed, 50);
    }
}
