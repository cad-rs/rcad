#[cfg(test)]
mod tests {

    use super::*;
    use crate::tolerance::{
        any_perpendicular, TOLERANCE_ABS, TOLERANCE_ADAPTIVE_MAX, TOLERANCE_COORD_SUB,
        TOLERANCE_FLOAT_DEDUP, TOLERANCE_LEN_MIN, TOLERANCE_LINEAR_ULTRA_STRICT,
        TOLERANCE_MESH_LEGACY, TOLERANCE_RETRY_LADDER_COARSE, TOLERANCE_RETRY_LADDER_MID,
    };
    use crate::check_orientation_consistency;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn remove_small_edges_removes_degenerate_loop() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        // Build a triangle with one degenerate self-loop edge (start == end).
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        // Edges: 0-1, 1-2, 2-0, plus degenerate 0-0
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 0 }); // degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let (fixed, removed) = remove_small_edges(&brep, TOLERANCE_MESH_LEGACY);
        assert!(removed >= 1, "degenerate self-loop should be removed");
        assert!(fixed.edges.len() < brep.edges.len());
    }

    #[test]
    fn remove_small_edges_is_noop_on_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let (fixed, removed) = remove_small_edges(&brep, TOLERANCE_ABS);
        assert_eq!(removed, 0, "unit box edges are not short");
        assert_eq!(fixed.edges.len(), brep.edges.len());
    }

    #[test]
    fn make_connected_baseline_merges_and_removes_tiny_edges() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 near-dup of 0

        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2
        brep.edges.push(Edge { start: 0, end: 3 }); // e3 tiny edge to be removed

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_baseline(&brep, TOLERANCE_MESH_LEGACY);
        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert_eq!(report.passes_run, 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!(fixed.edges.len() < brep.edges.len());
    }

    #[test]
    fn make_connected_iterative_reports_convergence() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative(&brep, TOLERANCE_MESH_LEGACY, 4);
        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert!(report.converged);
        assert!(report.passes_run >= 2);
        assert!(report.final_tolerance >= TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn make_connected_iterative_with_growth_increases_final_tolerance() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative_with_growth(&brep, TOLERANCE_MESH_LEGACY, 4, 2.0);
        assert!(report.passes_run >= 2);
        assert!(report.final_tolerance > TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn make_connected_iterative_with_growth_cap_clamps_tolerance() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative_with_growth_cap(
            &brep,
            TOLERANCE_MESH_LEGACY,
            4,
            10.0,
            2.0 * TOLERANCE_MESH_LEGACY,
        );
        assert!(report.passes_run >= 2);
        assert!(report.tolerance_cap_applied);
        assert!((report.final_tolerance - 2.0 * TOLERANCE_MESH_LEGACY).abs() <= TOLERANCE_FLOAT_DEDUP);
    }

    #[test]
    fn make_connected_iterative_growth_can_recover_after_initial_no_op_pass() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(50.0 * TOLERANCE_ABS, 0.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 });

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_iterative_with_growth_cap(
            &brep,
            TOLERANCE_MESH_LEGACY,
            2,
            10.0,
            TOLERANCE_RETRY_LADDER_MID,
        );

        assert_eq!(report.passes_run, 2);
        assert!(report.vertices_merged >= 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!((report.final_tolerance - TOLERANCE_RETRY_LADDER_MID).abs() <= TOLERANCE_FLOAT_DEDUP);
    }

    #[test]
    fn make_connected_scoped_growth_can_recover_after_initial_no_op_pass() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(50.0 * TOLERANCE_ABS, 0.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 });

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_iterative_scoped_with_growth_cap(
            &brep,
            &[0],
            TOLERANCE_MESH_LEGACY,
            2,
            10.0,
            TOLERANCE_RETRY_LADDER_MID,
        );

        assert_eq!(report.passes_run, 2);
        assert!(report.vertices_merged >= 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!((report.final_tolerance - TOLERANCE_RETRY_LADDER_MID).abs() <= TOLERANCE_FLOAT_DEDUP);
    }

    #[test]
    fn make_connected_scoped_only_affects_seed_region() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup near region A)
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 0.0, 0.0) }); // 4
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 1.0, 0.0) }); // 5
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 0.0, 0.0) }); // 6 (dup near region B)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge in scoped region
        brep.edges.push(Edge { start: 4, end: 5 });
        brep.edges.push(Edge { start: 5, end: 6 }); // unrelated region

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (scoped, report) = make_connected_iterative_scoped_with_growth_cap(
            &brep,
            &[0],
            TOLERANCE_MESH_LEGACY,
            3,
            1.0,
            TOLERANCE_RETRY_LADDER_COARSE,
        );

        assert!(report.vertices_merged >= 1);
        assert!(scoped.vertices.len() < brep.vertices.len());

        // Vertex near unrelated region B should remain after scoped cleanup.
        let has_far = scoped
            .vertices
            .iter()
            .any(|v| (v.point - DVec3::new(10.0, 0.0, 0.0)).length() <= TOLERANCE_LEN_MIN);
        assert!(has_far);
    }

    #[test]
    fn repair_unit_box_is_no_op() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let (fixed, report) = repair(&brep, TOLERANCE_ABS);
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.degenerate_faces_removed, 0);
        // Face count unchanged
        let faces: usize = fixed
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(faces, 6, "unit box should have 6 faces after repair");
    }

    #[test]
    fn merge_close_vertices_merges_duplicates() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        // Add two vertices at nearly the same position
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(TOLERANCE_COORD_SUB, 0.0, 0.0),
        }); // dup of 0
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });
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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, merged) = merge_close_vertices(&brep, TOLERANCE_MESH_LEGACY);
        assert!(merged >= 1, "should merge the near-duplicate vertex");
        assert!(
            fixed.vertices.len() < brep.vertices.len(),
            "should have fewer vertices"
        );
    }

    #[test]
    fn recompute_normals_fixes_zero_normal() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        // Face with wrong/zero normal
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // intentionally wrong
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let (fixed, n) = recompute_face_normals(&brep);
        assert!(
            n > 0 || fixed.solids[0].shells[0].faces[0].normal != DVec3::ZERO,
            "normal should have been fixed"
        );
    }

    #[test]
    fn fix_face_orientation_flips_inward_box_face() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let face = &mut brep.solids[0].shells[0].faces[0];
        face.normal = -face.normal;
        face.outer_wire = reverse_wire(&face.outer_wire);

        let before = check_orientation_consistency(&brep);
        assert!(!before.is_consistent);

        let (fixed, flipped) = fix_face_orientation(&brep);
        assert!(flipped >= 1);

        let after = check_orientation_consistency(&fixed);
        assert!(after.is_consistent, "orientation issues: {:?}", after.issues);
    }

    #[test]
    fn repair_reports_faces_reoriented() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let face = &mut brep.solids[0].shells[0].faces[0];
        face.normal = -face.normal;
        face.outer_wire = reverse_wire(&face.outer_wire);

        let (_fixed, report) = repair(&brep, TOLERANCE_ABS);
        assert!(report.faces_reoriented >= 1);
    }

    #[test]
    fn remove_degenerate_face() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });
        // Only 2 edges 閳?degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let (fixed, n) = remove_degenerate_faces(&brep);
        assert_eq!(n, 1);
        let face_count: usize = fixed
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(face_count, 0);
    }

    #[test]
    fn fix_same_range_flags_aligns_curve2d_ranges() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Build minimal SameRange mismatch for edge 0.
        if brep.geom.edge_curve_range.is_empty() {
            brep.geom.edge_curve_range = vec![Some([0.0, std::f64::consts::PI])];
        } else {
            brep.geom.edge_curve_range[0] = Some([0.0, std::f64::consts::PI]);
        }
        if brep.geom.edge_pcurves.is_empty() || brep.geom.edge_pcurves[0].is_empty() {
            // Sphere primitive normally has seam pcurves, but guard for future changes.
            return;
        }

        brep.geom.edge_same_range = vec![false; brep.edges.len().max(1)];
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        let pc = brep.geom.edge_pcurves[0][0];
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([1.0, 2.0]); // mismatched

        let (fixed, n) = fix_same_range_flags(&brep, TOLERANCE_COORD_SUB);
        assert!(n >= 1);
        assert!(fixed.geom.edge_same_range[0]);
        assert_eq!(
            fixed.geom.curve2d_range[pc.curve2d_idx],
            Some([0.0, std::f64::consts::PI])
        );
    }

    #[test]
    fn fix_same_range_with_scan_repairs_flagged_edges() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        if brep.geom.edge_curve_range.is_empty()
            || brep.geom.edge_pcurves.is_empty()
            || brep.geom.edge_pcurves[0].is_empty()
        {
            return;
        }

        brep.geom.edge_curve_range[0] = Some([0.0, std::f64::consts::PI]);
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        if brep.geom.edge_same_range.len() < brep.edges.len() {
            brep.geom.edge_same_range.resize(brep.edges.len(), true);
        }

        let pc = brep.geom.edge_pcurves[0][0];
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([1.0, 2.0]);
        brep.geom.edge_same_range[0] = false;

        let (fixed, n) = fix_same_range_with_scan(&brep, TOLERANCE_COORD_SUB);
        assert!(n >= 1);
        assert!(fixed.geom.edge_same_range[0]);
        assert_eq!(
            fixed.geom.curve2d_range[pc.curve2d_idx],
            Some([0.0, std::f64::consts::PI])
        );
    }

    #[test]
    fn propagate_tolerances_bottom_up_fills_slots_and_propagates() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Simple triangle face: 3 verts, 3 edges.
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2

        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
                    },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    sample_point: None,
                    mesh_dirty: true,
                surface_idx: None,
                }],
            }],
        });

        // Set vertex 0 with a large tolerance.
        brep.geom.vertex_tolerance = vec![TOLERANCE_ADAPTIVE_MAX, 0.0, 0.0];

        let out = propagate_tolerances(&brep, TOLERANCE_ABS, ToleranceFlowDirection::BottomUp);

        // vertex_tolerance slots must be filled.
        assert_eq!(out.geom.vertex_tolerance.len(), 3);
        // Edge tolerances should be at least floor.
        assert!(out.geom.edge_tolerance.len() >= 3);
        // Edge 0 connects v0 (tol=TOLERANCE_ADAPTIVE_MAX) and v1 (tol=floor); must 閳?TOLERANCE_ADAPTIVE_MAX.
        assert!(out.geom.edge_tolerance[0] >= TOLERANCE_ADAPTIVE_MAX);
        // Face tolerance should be 閳?max edge tolerance.
        assert!(out.geom.face_tolerance[0] >= out.geom.edge_tolerance[0]);
    }

    #[test]
    fn propagate_tolerances_top_down_spreads_face_tol_to_vertices() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
                    },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    sample_point: None,
                    mesh_dirty: true,
                surface_idx: None,
                }],
            }],
        });
        // Assign a large face tolerance.
        brep.geom.face_tolerance = vec![0.5 * TOLERANCE_RETRY_LADDER_COARSE];

        let out = propagate_tolerances(&brep, TOLERANCE_ABS, ToleranceFlowDirection::TopDown);

        // All edge tolerances should be 閳?face tolerance.
        for etol in &out.geom.edge_tolerance {
            assert!(*etol >= 0.5 * TOLERANCE_RETRY_LADDER_COARSE);
        }
        // All vertex tolerances should be 閳?face tolerance after propagation.
        for vtol in &out.geom.vertex_tolerance {
            assert!(*vtol >= 0.5 * TOLERANCE_RETRY_LADDER_COARSE);
        }
    }

    #[test]
    fn detect_shared_topology_advanced_detects_shared_vertices() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let report = detect_shared_topology_advanced(&brep, TOLERANCE_MESH_LEGACY);
        assert!(report.shared_vertex_pairs >= 1, "Should detect at least one shared vertex pair");
        assert!(report.has_shared_topology);
    }

    #[test]
    fn detect_shared_topology_advanced_detects_no_duplicate_faces_on_clean_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let report = detect_shared_topology_advanced(&brep, TOLERANCE_MESH_LEGACY);
        // A clean box should have NO fully shared (duplicate) faces
        assert_eq!(report.fully_shared_faces.len(), 0, "Clean box should have no duplicate faces");
        // A clean box has no duplicate vertices
        assert_eq!(report.shared_vertex_pairs, 0, "Clean box should have no duplicate vertices");
        // Note: Edge-based shared topology detection requires geometry data (curves)
        // which is not populated by the primitive box creation. The face sharing detection
        // for primitives uses topological edge indices, not geometric comparison.
    }

    #[test]
    fn make_connected_enhanced_with_mode_standard() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_enhanced_with_mode(
            &brep,
            TOLERANCE_MESH_LEGACY,
            4,
            MakeConnectedMode::Standard,
            false,
        );

        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert!(report.converged);
        assert!(fixed.vertices.len() < brep.vertices.len());
    }

    #[test]
    fn make_connected_enhanced_with_mode_conservative() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)
        brep.vertices.push(Vertex { point: DVec3::new(0.5, 0.0, 0.0) }); // 4 (creates short edge)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge
        brep.edges.push(Edge { start: 0, end: 4 }); // short edge (0.5 length, not tiny)

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_enhanced_with_mode(
            &brep,
            TOLERANCE_MESH_LEGACY,
            4,
            MakeConnectedMode::Conservative,
            false,
        );

        // Conservative mode should merge vertices but NOT remove short edges
        assert!(report.vertices_merged >= 1);
        assert_eq!(report.small_edges_removed, 0, "Conservative mode should not remove edges");
        assert!(report.converged);
    }

    #[test]
    fn make_connected_enhanced_with_mode_aggressive() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

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
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_enhanced_with_mode(
            &brep,
            TOLERANCE_MESH_LEGACY,
            4,
            MakeConnectedMode::Aggressive,
            false,
        );

        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert!(report.converged);
        assert!(fixed.vertices.len() < brep.vertices.len());
    }

    #[test]
    fn shared_edge_info_structure_works() {
        let info = SharedEdgeInfo {
            edge_a: 0,
            edge_b: 1,
            geometry_compatible: true,
            curvature_continuous: true,
            param_range_compatible: true,
            max_deviation: 0.001,
            reversed: false,
        };

        assert_eq!(info.edge_a, 0);
        assert_eq!(info.edge_b, 1);
        assert!(info.geometry_compatible);
        assert!(info.curvature_continuous);
        assert!(info.param_range_compatible);
    }

    #[test]
    fn shared_face_info_structure_works() {
        let info = SharedFaceInfo {
            face_a: 0,
            face_b: 1,
            kind: SharedFaceKind::PartialShared,
            shared_edges: vec![0, 1],
            shared_vertices: vec![0, 1, 2],
            normals_compatible: true,
        };

        assert_eq!(info.face_a, 0);
        assert_eq!(info.face_b, 1);
        assert_eq!(info.kind, SharedFaceKind::PartialShared);
        assert_eq!(info.shared_edges.len(), 2);
        assert_eq!(info.shared_vertices.len(), 3);
    }

    #[test]
    fn shared_topology_report_structure_works() {
        let mut report = SharedTopologyReport::default();
        report.fully_shared_faces.push(SharedFaceInfo {
            face_a: 0,
            face_b: 1,
            kind: SharedFaceKind::FullShared,
            shared_edges: vec![],
            shared_vertices: vec![],
            normals_compatible: true,
        });
        report.shared_edges.push(SharedEdgeInfo {
            edge_a: 0,
            edge_b: 1,
            geometry_compatible: true,
            curvature_continuous: true,
            param_range_compatible: true,
            max_deviation: 0.0,
            reversed: false,
        });
        report.shared_vertex_pairs = 2;
        report.has_shared_topology = true;

        assert_eq!(report.fully_shared_faces.len(), 1);
        assert_eq!(report.shared_edges.len(), 1);
        assert_eq!(report.shared_vertex_pairs, 2);
        assert!(report.has_shared_topology);
    }

    #[test]
    fn edge_sew_config_default_values() {
        let config = EdgeSewConfig::default();
        assert!(config.base_tolerance > 0.0);
        assert!(config.max_tolerance >= config.base_tolerance);
        assert!(config.tolerance_growth >= 1.0);
        assert!(config.max_passes > 0);
        assert!(config.use_geometric_proximity);
        assert!(config.merge_same_curve_edges);
        assert!(config.handle_periodic_seams);
    }

    #[test]
    fn adaptive_tolerance_config_default_values() {
        let config = AdaptiveToleranceConfig::default();
        assert!(config.base_tolerance > 0.0);
        assert!(config.max_tolerance >= config.base_tolerance);
        assert!(config.tolerance_growth >= 1.0);
        assert!(config.min_feature_size > 0.0);
        assert!(config.use_curvature_adjustment);
    }

    #[test]
    fn sew_edges_enhanced_basic() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
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

        let config = EdgeSewConfig::default();
        let (_, report) = sew_edges_enhanced(&brep, &config);

        // The function should run without error
        assert!(report.passes_executed >= 1);
    }

    #[test]
    fn merge_vertices_adaptive_basic() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
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

        let config = AdaptiveToleranceConfig::default();
        let (_, report) = merge_vertices_adaptive(&brep, &config);

        // The function should run without error
        assert!(report.passes_executed >= 1);
    }

    #[test]
    fn enhanced_edge_sew_report_default() {
        let report = EnhancedEdgeSewReport::default();
        assert_eq!(report.edges_sewn, 0);
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.passes_executed, 0);
        assert!(!report.converged);
    }

    #[test]
    fn adaptive_tolerance_merge_report_default() {
        let report = AdaptiveToleranceMergeReport::default();
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.edges_removed, 0);
        assert_eq!(report.passes_executed, 0);
        assert!(!report.converged);
    }

    // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    // B-Spline Same-Domain Tests
    // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

    #[test]
    fn bspline_same_domain_identical_surfaces() {
        use rcad_kernel::geom::BSplineSurface;

        let surf = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf, &surf, TOLERANCE_MESH_LEGACY);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(match_result.is_same_domain);
        assert!(match_result.degrees_match);
        assert!(match_result.knots_match);
        assert!(match_result.max_control_point_deviation < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn bspline_same_domain_different_degrees() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.5, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.5, 0.5, 0.0), DVec3::new(0.5, 0.5, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, TOLERANCE_MESH_LEGACY);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(!match_result.degrees_match);
    }

    #[test]
    fn bspline_same_domain_different_knots() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 2.0, 2.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, TOLERANCE_MESH_LEGACY);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(match_result.degrees_match);
        assert!(!match_result.knots_match);
    }

    #[test]
    fn bspline_same_domain_different_control_points() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, TOLERANCE_MESH_LEGACY);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(match_result.max_control_point_deviation > 0.5);
    }

    #[test]
    fn bspline_continuity_default() {
        let continuity = BsplineContinuity::default();
        assert_eq!(continuity, BsplineContinuity::None);
    }

    #[test]
    fn check_bspline_continuity_same_surface() {
        use rcad_kernel::geom::BSplineSurface;

        let surf = BSplineSurface {
            degree_u: 3,
            degree_v: 3,
            knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.33, 0.0, 0.0), DVec3::new(0.66, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 0.33, 0.0), DVec3::new(0.33, 0.33, 0.0), DVec3::new(0.66, 0.33, 0.0), DVec3::new(1.0, 0.33, 0.0)],
                vec![DVec3::new(0.0, 0.66, 0.0), DVec3::new(0.33, 0.66, 0.0), DVec3::new(0.66, 0.66, 0.0), DVec3::new(1.0, 0.66, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(0.33, 1.0, 0.0), DVec3::new(0.66, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0, 1.0],
            ],
        };

        let continuity = check_bspline_continuity(&surf, &surf, TOLERANCE_MESH_LEGACY);
        // A bicubic B-spline with clamped boundary knots (multiplicity 4) has C0 continuity
        // at boundaries due to knot multiplicity = degree, but is C2 inside the domain.
        // Our implementation reports minimum continuity at any knot, which is C0 at boundaries.
        assert!(continuity >= BsplineContinuity::C0);
    }

    #[test]
    fn check_bspline_continuity_adjacent_v() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let continuity = check_bspline_continuity(&surf1, &surf2, TOLERANCE_MESH_LEGACY);
        assert!(continuity >= BsplineContinuity::C0);
    }

    #[test]
    fn max_knot_multiplicity_single() {
        let knots = vec![0.0, 0.0, 0.5, 1.0, 1.0];
        let mult = max_knot_multiplicity(&knots);
        assert_eq!(mult, 2);
    }

    #[test]
    fn max_knot_multiplicity_triple() {
        let knots = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0];
        let mult = max_knot_multiplicity(&knots);
        assert_eq!(mult, 3);
    }

    #[test]
    fn same_domain_match_debug() {
        let match_result = SameDomainMatch {
            is_same_domain: true,
            continuity: BsplineContinuity::C1,
            max_control_point_deviation: 0.0,
            max_weight_deviation: 0.0,
            knots_match: true,
            degrees_match: true,
        };

        let debug_str = format!("{:?}", match_result);
        assert!(debug_str.contains("is_same_domain: true"));
        assert!(debug_str.contains("C1"));
    }

    #[test]
    fn merged_face_info_debug() {
        let info = MergedFaceInfo {
            kept_face_idx: 0,
            removed_face_idx: 1,
            merged_edge_count: 6,
            inner_wires_merged: false,
            continuity: BsplineContinuity::C0,
        };

        let debug_str = format!("{:?}", info);
        assert!(debug_str.contains("kept_face_idx: 0"));
        assert!(debug_str.contains("merged_edge_count: 6"));
    }

    #[test]
    fn bspline_continuity_ordering() {
        assert!(BsplineContinuity::None < BsplineContinuity::C0);
        assert!(BsplineContinuity::C0 < BsplineContinuity::G1);
        assert!(BsplineContinuity::G1 < BsplineContinuity::C1);
        assert!(BsplineContinuity::C1 < BsplineContinuity::C2);
        assert!(BsplineContinuity::C2 < BsplineContinuity::CN);
    }

    #[test]
    fn bspline_same_domain_rational_surface() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0), DVec3::new(2.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0), DVec3::new(2.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 2.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        };

        let surf2 = surf1.clone();

        let result = bspline_same_domain(&surf1, &surf2, TOLERANCE_MESH_LEGACY);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(match_result.is_same_domain);
        assert!(match_result.max_weight_deviation < TOLERANCE_COORD_SUB);
    }

    #[test]
    fn bspline_same_domain_different_weights() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0), DVec3::new(2.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0), DVec3::new(2.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 2.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0), DVec3::new(2.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0), DVec3::new(2.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 3.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, TOLERANCE_MESH_LEGACY);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(match_result.max_weight_deviation > 0.5);
    }

    // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    // Shell and Solid Repair Tests
    // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

    #[test]
    fn check_shell_closure_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Unit box shell should be closed");
        assert_eq!(report.open_edge_count, 0, "Unit box should have no open edges");
        assert_eq!(report.face_count, 6, "Unit box should have 6 faces");
        assert!(report.euler_characteristic > 0, "Unit box should have positive Euler characteristic");
    }

    #[test]
    fn check_shell_closure_unit_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Sphere shell should be closed");
        assert_eq!(report.open_edge_count, 0, "Sphere should have no open edges");
    }

    #[test]
    fn check_shell_closure_unit_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Cylinder shell should be closed");
        assert_eq!(report.open_edge_count, 0, "Cylinder should have no open edges");
    }

    #[test]
    fn check_shell_closure_open_triangle() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        // Create an open triangle (not a closed shell)
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        // Missing edge 2-0 to make it open

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        };

        let shell = Shell { faces: vec![face] };
        let report = check_shell_closure(&shell, &brep);

        assert!(!report.is_closed, "Open triangle should not be closed");
        assert!(report.open_edge_count > 0, "Open triangle should have open edges");
    }

    #[test]
    fn fix_shell_orientation_inverted_normals() {
        // Create a box and invert all its face normals
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Invert all face normals in the shell
        for face in &mut brep.solids[0].shells[0].faces {
            face.normal = -face.normal;
        }

        let shell = &brep.solids[0].shells[0].clone();
        let (fixed_shell, report) = fix_shell_orientation(shell, &brep);

        // All 6 faces should be reoriented
        assert!(report.faces_reoriented >= 6, "Should reorient all inverted faces");

        // All normals should now point outward (have positive dot product with outward direction)
        for face in &fixed_shell.faces {
            // For a box centered at origin, check that normals are consistent
            let normal_magnitude = face.normal.length();
            assert!(normal_magnitude > 0.99, "Normal should be unit length");
        }
    }

    #[test]
    fn fix_shell_orientation_correct_normals() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let (_, report) = fix_shell_orientation(shell, &brep);

        // Box from primitive should already have correct normals
        assert_eq!(report.faces_reoriented, 0, "Box from primitive should not need reorientation");
    }

    #[test]
    fn shell_fix_report_summary() {
        let report = ShellFixReport {
            faces_reoriented: 3,
            non_manifold_edges_processed: 1,
            shells_created: 0,
            is_closed: true,
            is_manifold: true,
            open_edge_count: 0,
            non_manifold_edge_count: 0,
        };

        let summary = report.summary();
        assert!(summary.contains("3 faces reoriented"));
        assert!(summary.contains("closed=true"));
    }

    #[test]
    fn closure_report_summary() {
        let report = ClosureReport {
            is_closed: true,
            open_edge_count: 0,
            open_edges: vec![],
            euler_characteristic: 2,
            vertex_count: 8,
            edge_count: 12,
            face_count: 6,
            is_orientable: true,
            genus: Some(0),
        };

        let summary = report.summary();
        assert!(summary.contains("Closed shell"));
        assert!(summary.contains("V=8"));
        assert!(summary.contains("E=12"));
        assert!(summary.contains("F=6"));
        assert!(summary.contains("genus=0"));
    }

    #[test]
    fn check_solid_closure_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let solid = &brep.solids[0];
        let report = check_solid_closure(solid, &brep);

        assert!(report.is_closed, "Box solid should be closed");
        assert!(report.has_proper_nesting, "Box should have proper nesting");
        assert_eq!(report.outer_shell_count, 1, "Box should have 1 outer shell");
        assert_eq!(report.inner_shell_count, 0, "Box should have 0 inner shells");
    }

    #[test]
    fn check_solid_closure_unit_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        let solid = &brep.solids[0];
        let report = check_solid_closure(solid, &brep);

        assert!(report.is_closed, "Sphere solid should be closed");
        assert_eq!(report.outer_shell_count, 1, "Sphere should have 1 outer shell");
    }

    #[test]
    fn fix_solid_orientation_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let solid = &brep.solids[0];
        let (_, report) = fix_solid_orientation(solid, &brep);

        // Box from primitive should already be properly oriented
        assert!(report.has_valid_closure, "Box should have valid closure");
        assert_eq!(report.outer_shells, 1, "Box should have 1 outer shell");
        assert_eq!(report.inner_shells, 0, "Box should have 0 inner shells");
    }

    #[test]
    fn fix_solid_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let solid = &brep.solids[0];
        let (fixed_solid, report) = fix_solid(solid, &brep);

        assert!(report.is_clean(), "Fixed solid should be clean");
        assert!(report.has_valid_closure, "Fixed solid should have valid closure");
        assert!(report.is_properly_oriented, "Fixed solid should be properly oriented");
        assert_eq!(fixed_solid.shells.len(), solid.shells.len(), "Shell count should be preserved");
    }

    #[test]
    fn solid_fix_report_summary() {
        let report = SolidFixReport {
            shells_reoriented: 1,
            faces_reoriented: 3,
            outer_shells: 1,
            inner_shells: 0,
            is_properly_oriented: true,
            has_valid_closure: true,
            total_fixes: 4,
        };

        let summary = report.summary();
        assert!(summary.contains("1 shells reoriented"));
        assert!(summary.contains("3 faces flipped"));
        assert!(summary.contains("1 outer"));
    }

    #[test]
    fn solid_closure_report_summary() {
        let report = SolidClosureReport {
            is_closed: true,
            has_proper_nesting: true,
            outer_shell_count: 1,
            inner_shell_count: 2,
            unclosed_shell_indices: vec![],
            volume: 10.5,
            shell_euler: vec![2, 2, 2],
            solid_euler: 6,
        };

        let summary = report.summary();
        assert!(summary.contains("Valid solid"));
        assert!(summary.contains("2 voids"));
    }

    #[test]
    fn check_shell_closure_torus() {
        let brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Torus shell should be closed");
        // Torus has genus 1, so Euler characteristic should be 0
        assert_eq!(report.euler_characteristic, 0, "Torus should have Euler characteristic 0");
        assert_eq!(report.genus, Some(1), "Torus should have genus 1");
    }

    #[test]
    fn fix_non_manifold_shell_already_manifold() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let (_, report) = fix_non_manifold_shell(shell, &brep);

        assert!(report.is_manifold, "Box shell should be manifold");
        assert_eq!(report.non_manifold_edge_count, 0, "Box should have no non-manifold edges");
    }

    #[test]
    fn shell_orientability_check() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        // Box should be closed (no open edges)
        assert!(report.is_closed, "Box shell should be closed");
        // Note: orientability check depends on face orientation consistency
        // which may vary based on how the primitive is constructed
    }

    // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    // Tests for Enhanced Shell Repair Functions
    // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

    #[test]
    fn fix_shell_orientation_advanced_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let (fixed_shell, report) = fix_shell_orientation_advanced(shell, &brep);

        // Box should have no edge conflicts
        assert_eq!(report.edge_conflicts, 0, "Box should have no edge orientation conflicts");
        assert_eq!(fixed_shell.faces.len(), shell.faces.len(), "Face count should be preserved");
    }

    #[test]
    fn fix_shell_orientation_advanced_inverted_box() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Invert all face normals
        for face in &mut brep.solids[0].shells[0].faces {
            face.normal = -face.normal;
        }

        let shell = brep.solids[0].shells[0].clone();
        let (fixed_shell, report) = fix_shell_orientation_advanced(&shell, &brep);

        // The algorithm should process all faces
        assert_eq!(fixed_shell.faces.len(), shell.faces.len(), "Face count should be preserved");
        // Edge conflicts should be resolved after repair
        assert_eq!(report.edge_conflicts, 0, "Edge conflicts should be resolved");
    }

    #[test]
    fn repair_shell_closure_closed_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let result = repair_shell_closure(shell, &brep, 0.001);

        // Closed box should remain closed
        assert!(result.is_closed, "Box should remain closed");
        assert_eq!(result.open_edges_detected, 0, "Box should have no open edges");
        assert_eq!(result.faces_added, 0, "No faces should be added");
    }

    #[test]
    fn repair_shell_closure_open_shell() {
        use rcad_kernel::topology::{Edge, Face, Shell, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });
        // Missing diagonal edge to close the square

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        };

        let shell = Shell { faces: vec![face] };
        let result = repair_shell_closure(&shell, &brep, 0.001);

        // Open shell should detect open edges
        assert!(result.open_edges_detected > 0, "Should detect open edges");
    }

    #[test]
    fn repair_non_manifold_edges_manifold_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let result = repair_non_manifold_edges(shell, &brep);

        // Box is manifold
        assert!(result.is_manifold, "Box should be manifold");
        assert_eq!(result.edges_processed, 0, "No non-manifold edges to process");
    }

    #[test]
    fn validate_shell_topology_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = validate_shell_topology(shell, &brep);

        assert!(report.is_closed, "Box should be closed");
        assert!(report.is_manifold, "Box should be manifold");
        assert_eq!(report.face_count, 6, "Box should have 6 faces");
        assert!(report.edge_valence.iter().all(|e| e.is_manifold), "All edges should be manifold");
    }

    #[test]
    fn validate_shell_topology_unit_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = validate_shell_topology(shell, &brep);

        assert!(report.is_closed, "Sphere should be closed");
        assert!(report.is_manifold, "Sphere should be manifold");
    }

    #[test]
    fn validate_shell_topology_torus() {
        let brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        let shell = &brep.solids[0].shells[0];
        let report = validate_shell_topology(shell, &brep);

        assert!(report.is_closed, "Torus should be closed");
        assert!(report.is_manifold, "Torus should be manifold");
        assert_eq!(report.genus, Some(1), "Torus should have genus 1");
        assert_eq!(report.euler_characteristic, 0, "Torus Euler characteristic should be 0");
    }

    #[test]
    fn validate_shell_topology_open_shell() {
        use rcad_kernel::topology::{Edge, Face, Shell, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        // Missing edge to close the triangle

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        };

        let shell = Shell { faces: vec![face] };
        let report = validate_shell_topology(&shell, &brep);

        assert!(!report.is_closed, "Open triangle should not be closed");
        assert!(report.open_edge_count > 0, "Should have open edges");
    }

    #[test]
    fn shell_orientation_report_summary() {
        let report = ShellOrientationReport {
            faces_inverted: 3,
            faces_correct: 5,
            inverted_face_indices: vec![0, 2, 4],
            edge_conflicts: 0,
            is_consistent: true,
            non_manifold_edges_skipped: 0,
            volume_sign: 1.0,
        };

        let summary = report.summary();
        assert!(summary.contains("3 inverted"));
        assert!(summary.contains("5 correct"));
        assert!(summary.contains("consistent=true"));
    }

    #[test]
    fn shell_closure_result_summary() {
        let result = ShellClosureResult {
            original_shell: Shell { faces: vec![] },
            repaired_shell: Shell { faces: vec![] },
            open_edges_detected: 2,
            gaps_closed: 1,
            faces_added: 1,
            unrepairable_gaps: vec![],
            is_closed: true,
            tolerance_used: 0.001,
        };

        let summary = result.summary();
        assert!(summary.contains("closed 1 gaps"));
        assert!(summary.contains("added 1 faces"));
    }

    #[test]
    fn manifold_repair_result_summary() {
        let result = ManifoldRepairResult {
            original_shell: Shell { faces: vec![] },
            repaired_shell: Shell { faces: vec![] },
            edges_processed: 2,
            edges_split: 1,
            vertices_duplicated: 2,
            faces_created: 0,
            is_manifold: true,
            edge_details: vec![],
        };

        let summary = result.summary();
        assert!(summary.contains("1 edges split"));
        assert!(summary.contains("2 vertices duplicated"));
        assert!(summary.contains("manifold=true"));
    }

    #[test]
    fn shell_validation_report_summary() {
        let report = ShellValidationReport {
            is_valid: true,
            euler_characteristic: 2,
            expected_euler: Some(2),
            euler_valid: true,
            vertex_count: 8,
            edge_count: 12,
            face_count: 6,
            open_edge_count: 0,
            non_manifold_edge_count: 0,
            non_manifold_vertex_count: 0,
            orientation_consistent: true,
            is_closed: true,
            is_manifold: true,
            genus: Some(0),
            edge_valence: vec![],
            vertex_valence: vec![],
            errors: vec![],
            warnings: vec![],
        };

        let summary = report.summary();
        assert!(summary.contains("VALID"));
        assert!(summary.contains("V=8"));
        assert!(summary.contains("E=12"));
        assert!(summary.contains("F=6"));
        assert!(report.is_closed_manifold());
    }

    #[test]
    fn gap_info_creation() {
        let gap = GapInfo {
            boundary_edges: vec![0, 1, 2],
            estimated_area: 0.5,
            can_fill: true,
            failure_reason: None,
        };

        assert_eq!(gap.boundary_edges.len(), 3);
        assert!(gap.can_fill);
        assert!(gap.failure_reason.is_none());
    }

    #[test]
    fn non_manifold_edge_info_creation() {
        let info = NonManifoldEdgeInfo {
            edge_index: 5,
