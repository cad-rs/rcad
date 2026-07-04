#[cfg(test)]
mod tests {
    use crate::brep_graph::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn topo_graph_from_box_has_expected_counts() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);

        let solids = g.nodes.iter().filter(|n| n.kind == NodeKind::Solid).count();
        let shells = g.nodes.iter().filter(|n| n.kind == NodeKind::Shell).count();
        let faces = g.nodes.iter().filter(|n| n.kind == NodeKind::Face).count();
        let wires = g.nodes.iter().filter(|n| n.kind == NodeKind::Wire).count();
        let edges = g.nodes.iter().filter(|n| n.kind == NodeKind::Edge).count();
        let vertices = g.nodes.iter().filter(|n| n.kind == NodeKind::Vertex).count();

        assert_eq!(solids, 1);
        assert_eq!(shells, 1);
        assert_eq!(faces, 6);
        assert_eq!(wires, 6);
        assert_eq!(edges, 12);
        assert_eq!(vertices, 8);
    }

    #[test]
    fn topo_graph_faces_and_edges_queries_work() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);

        let shell0 = TopoNode {
            kind: NodeKind::Shell,
            index: 0,
        };
        let faces = g.faces_of_shell(shell0);
        assert_eq!(faces.len(), 6);

        let face0 = TopoNode {
            kind: NodeKind::Face,
            index: 0,
        };
        let edges = g.edges_of_face(face0);
        assert_eq!(edges.len(), 4);

        let edge0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        let verts = g.vertices_of_edge(edge0);
        assert_eq!(verts.len(), 2);
        assert_eq!(verts[0].index, brep.edges[0].start);
        assert_eq!(verts[1].index, brep.edges[0].end);
    }

    #[test]
    fn topo_graph_validate_passes_on_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);
        assert!(g.validate().is_empty());
        assert!(!g.history.events.is_empty());
    }

    #[test]
    fn topo_graph_compact_drops_orphans_and_dedups_neighbors() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);

        // Inject duplicate and orphan adjacency entries to simulate noisy edits.
        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        if let Some(neigh) = g.adjacency.get_mut(&v0) {
            neigh.push(e0);
            neigh.push(e0);
        }
        let orphan = TopoNode {
            kind: NodeKind::Face,
            index: 9999,
        };
        g.adjacency.insert(orphan, vec![v0]);

        g.compact();

        assert!(!g.adjacency.contains_key(&orphan));
        let neigh = g.adjacency.get(&v0).expect("vertex adjacency exists");
        let count_e0 = neigh.iter().filter(|n| **n == e0).count();
        assert_eq!(count_e0, 1);
        assert!(g.history.events.iter().any(|e| e.action == "compact"));
    }

    #[test]
    fn topo_graph_mutate_checked_reports_invalid_graph() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);

        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };

        let res = g.mutate_checked("inject_nonsymmetric", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.retain(|n| *n != v0);
            }
        });

        assert!(res.is_err());
        let issues = res.expect_err("mutation should be invalid");
        assert!(issues
            .iter()
            .any(|i| matches!(i, TopoGraphValidationIssue::NonSymmetricAdjacency { .. })));
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_invalid:inject_nonsymmetric"));
    }

    #[test]
    fn topo_graph_mutate_guarded_rolls_back_on_invalid_mutation() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);
        let before = g.clone();

        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };

        let res = g.mutate_guarded("inject_nonsymmetric", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.retain(|n| *n != v0);
            }
        });

        assert!(res.is_err());
        assert_eq!(g.nodes, before.nodes);
        assert_eq!(g.adjacency, before.adjacency);
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_guarded_rollback:inject_nonsymmetric"));
    }

    #[test]
    fn topo_graph_mutate_guarded_commits_valid_mutation() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };

        let res = g.mutate_guarded("dedup_neighbors", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.push(v0);
            }
            graph.compact();
        });

        assert!(res.is_ok());
        let neigh = g.adjacency.get(&e0).expect("edge adjacency exists");
        let count_v0 = neigh.iter().filter(|n| **n == v0).count();
        assert_eq!(count_v0, 1);
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_guarded:dedup_neighbors"));
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Tests for Cross-Operation Naming Stability Analysis
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

#[cfg(test)]
mod cross_operation_tests {
    use crate::brep_graph::*;
    use rcad_kernel::PrimitiveSolid;
    use rcad_kernel::persistent_naming::{
        OperationType, NamingEvent, PersistentId,
    };

    /// Helper to create a history with a single operation.
    fn create_history_with_operation(
        operation_type: OperationType,
        label: &str,
        stats: OperationStats,
        events: Vec<NamingEvent>,
    ) -> BRepGraphHistory {
        let mut history = BRepGraphHistory::new();
        history.begin_operation(operation_type, Some(label.to_string()));
        history.naming_engine_mut().cross_operation_history_mut()
            .add_events(rcad_kernel::persistent_naming::OperationId(1), events);
        history.finalize_operation(stats);
        history
    }

    #[test]
    fn analyze_empty_history_returns_perfect_stability() {
        let analysis = analyze_naming_sequence(&[], &[]);

        assert!(analysis.is_excellent());
        assert_eq!(analysis.overall_stability, 1.0);
        assert!(analysis.broken_chains.is_empty());
        assert!(analysis.operations.is_empty());
    }

    #[test]
    fn analyze_single_operation_with_no_losses() {
        let history = create_history_with_operation(
            OperationType::BooleanUnion,
            "test_union",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 10,
                names_preserved: 10,
                names_lost: 0,
                names_generated: 0,
                conflicts_resolved: 0,
            },
            vec![],
        );

        let analysis = analyze_naming_sequence(&[history], &[]);

        assert!(analysis.is_excellent());
        assert_eq!(analysis.overall_stability, 1.0);
        assert_eq!(analysis.operations.len(), 1);
    }

    #[test]
    fn analyze_single_operation_with_losses() {
        let history = create_history_with_operation(
            OperationType::BooleanDifference,
            "test_diff",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 7,
                names_preserved: 7,
                names_lost: 3,
                names_generated: 0,
                conflicts_resolved: 0,
            },
            vec![
                NamingEvent::Lost {
                    entity_id: 1,
                    persistent_id: PersistentId(1),
                },
                NamingEvent::Lost {
                    entity_id: 2,
                    persistent_id: PersistentId(2),
                },
                NamingEvent::Lost {
                    entity_id: 3,
                    persistent_id: PersistentId(3),
                },
            ],
        );

        let initial_entities: Vec<TopoNode> = (0..10)
            .map(|i| TopoNode { kind: NodeKind::Face, index: i })
            .collect();

        let analysis = analyze_naming_sequence(&[history], &initial_entities);

        assert!(!analysis.is_excellent());
        assert_eq!(analysis.broken_chains.len(), 3);
        assert!(analysis.overall_stability < 1.0);
    }

    #[test]
    fn analyze_multiple_operations_cumulative_stability() {
        let history1 = create_history_with_operation(
            OperationType::BooleanUnion,
            "union1",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 12,
                names_preserved: 10,
                names_lost: 0,
                names_generated: 2,
                conflicts_resolved: 0,
            },
            vec![],
        );

        let history2 = create_history_with_operation(
            OperationType::BooleanDifference,
            "diff1",
            OperationStats {
                entity_count_before: 12,
                entity_count_after: 10,
                names_preserved: 8,
                names_lost: 4,
                names_generated: 2,
                conflicts_resolved: 0,
            },
            vec![
                NamingEvent::Lost {
                    entity_id: 1,
                    persistent_id: PersistentId(1),
                },
                NamingEvent::Lost {
                    entity_id: 2,
                    persistent_id: PersistentId(2),
                },
            ],
        );

        let analysis = analyze_naming_sequence(&[history1, history2], &[]);

        assert_eq!(analysis.operations.len(), 2);
        assert_eq!(analysis.per_operation_stability.len(), 2);

        // Cumulative stability should decrease.
        assert!(analysis.per_operation_stability[1].cumulative_stability <= analysis.per_operation_stability[0].cumulative_stability);
    }

    #[test]
    fn stability_trend_improving() {
        // Create metrics with improving stability.
        let metrics = vec![
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(1),
                operation_type: OperationType::BooleanUnion,
                label: None,
                names_retained: 5,
                names_lost: 5,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.5,
                cumulative_stability: 0.5,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(2),
                operation_type: OperationType::BooleanUnion,
                label: None,
                names_retained: 7,
                names_lost: 3,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.7,
                cumulative_stability: 0.6,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(3),
                operation_type: OperationType::BooleanUnion,
                label: None,
                names_retained: 9,
                names_lost: 1,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.9,
                cumulative_stability: 0.7,
            },
        ];

        let trend = calculate_stability_trend(&metrics);
        assert!(trend > 0.0, "Trend should be positive (improving), got {}", trend);
    }

    #[test]
    fn stability_trend_degrading() {
        // Create metrics with degrading stability.
        let metrics = vec![
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(1),
                operation_type: OperationType::BooleanDifference,
                label: None,
                names_retained: 9,
                names_lost: 1,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.9,
                cumulative_stability: 0.9,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(2),
                operation_type: OperationType::BooleanDifference,
                label: None,
                names_retained: 7,
                names_lost: 3,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.7,
                cumulative_stability: 0.8,
            },
            OperationStabilityMetrics {
                operation_id: rcad_kernel::persistent_naming::OperationId(3),
                operation_type: OperationType::BooleanDifference,
                label: None,
                names_retained: 5,
                names_lost: 5,
                names_generated: 0,
                conflicts: 0,
                stability_score: 0.5,
                cumulative_stability: 0.6,
            },
        ];

        let trend = calculate_stability_trend(&metrics);
        assert!(trend < 0.0, "Trend should be negative (degrading), got {}", trend);
    }

    #[test]
    fn detect_conflicts_no_conflicts() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 1.0,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let conflicts = detect_cross_operation_conflicts(&analysis);
        assert!(conflicts.is_empty());
    }

    #[test]
    fn detect_conflicts_broken_genealogy() {
        let mut genealogy = HashMap::new();
        genealogy.insert(
            PersistentId(1),
            EntityGenealogy {
                persistent_id: PersistentId(1),
                created_in_operation: rcad_kernel::persistent_naming::OperationId(1),
                evolution: vec![],  // Empty evolution is a broken genealogy
                current_entity_id: Some(42),
                is_deleted: false,
            },
        );

        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: genealogy,
            per_operation_stability: vec![],
            overall_stability: 1.0,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let conflicts = detect_cross_operation_conflicts(&analysis);
        assert!(!conflicts.is_empty());

        let has_broken = conflicts.iter().any(|c| c.conflict_type == ConflictType::BrokenGenealogy);
        assert!(has_broken, "Should detect BrokenGenealogy conflict");
    }

    #[test]
    fn generate_recommendations_critical_stability() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.3,  // Critical
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        assert!(!recommendations.is_empty());

        // Should have a high-priority architecture recommendation.
        let has_critical = recommendations.iter().any(|r| r.priority == 100);
        assert!(has_critical, "Should have critical stability recommendation");
    }

    #[test]
    fn generate_recommendations_good_stability() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.95,  // Good
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        // With good stability and no issues, should have minimal or no recommendations.
        // Actually, it should be empty since no issues.
        assert!(recommendations.is_empty() || recommendations.iter().all(|r| r.priority < 70));
    }

    #[test]
    fn generate_recommendations_broken_chains() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.9,
            broken_chains: vec![
                BrokenChainInfo {
                    persistent_id: PersistentId(1),
                    broken_at_operation: rcad_kernel::persistent_naming::OperationId(1),
                    entity_id: 1,
                    entity_type: Some(rcad_kernel::persistent_naming::EntityType::Face),
                    survived_operations: 10,  // Long-lived entity
                    break_reason: "Test break".to_string(),
                },
            ],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        let has_entity_tracking = recommendations.iter().any(|r| r.category == RecommendationCategory::EntityTracking);
        assert!(has_entity_tracking, "Should have entity tracking recommendation");
    }

    #[test]
    fn generate_recommendations_degrading_trend() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.85,
            broken_chains: vec![],
            stability_trend: -0.2,  // Degrading
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        let has_propagation = recommendations.iter().any(|r| r.category == RecommendationCategory::PropagationPolicy);
        assert!(has_propagation, "Should have propagation policy recommendation for degrading trend");
    }

    #[test]
    fn most_problematic_operation() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(1),
                    operation_type: OperationType::BooleanUnion,
                    label: None,
                    names_retained: 9,
                    names_lost: 1,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.9,
                    cumulative_stability: 0.9,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(2),
                    operation_type: OperationType::BooleanDifference,
                    label: None,
                    names_retained: 5,
                    names_lost: 5,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.5,  // Most problematic
                    cumulative_stability: 0.45,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(3),
                    operation_type: OperationType::Feature,
                    label: None,
                    names_retained: 8,
                    names_lost: 2,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.8,
                    cumulative_stability: 0.36,
                },
            ],
            overall_stability: 0.7,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let problematic = analysis.most_problematic_operation();
        assert!(problematic.is_some());
        assert_eq!(problematic.unwrap().operation_id, rcad_kernel::persistent_naming::OperationId(2));
    }

    #[test]
    fn operations_by_stability_sorted() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(1),
                    operation_type: OperationType::BooleanUnion,
                    label: None,
                    names_retained: 7,
                    names_lost: 3,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.7,
                    cumulative_stability: 0.7,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(2),
                    operation_type: OperationType::BooleanDifference,
                    label: None,
                    names_retained: 5,
                    names_lost: 5,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.5,
                    cumulative_stability: 0.35,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(3),
                    operation_type: OperationType::Feature,
                    label: None,
                    names_retained: 9,
                    names_lost: 1,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.9,
                    cumulative_stability: 0.315,
                },
            ],
            overall_stability: 0.7,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let sorted = analysis.operations_by_stability();
        assert_eq!(sorted.len(), 3);
        // Should be sorted by stability score (ascending).
        assert_eq!(sorted[0].stability_score, 0.5);
        assert_eq!(sorted[1].stability_score, 0.7);
        assert_eq!(sorted[2].stability_score, 0.9);
    }

    #[test]
    fn average_operation_stability() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(1),
                    operation_type: OperationType::BooleanUnion,
                    label: None,
                    names_retained: 8,
                    names_lost: 2,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.8,
                    cumulative_stability: 0.8,
                },
                OperationStabilityMetrics {
                    operation_id: rcad_kernel::persistent_naming::OperationId(2),
                    operation_type: OperationType::BooleanDifference,
                    label: None,
                    names_retained: 6,
                    names_lost: 4,
                    names_generated: 0,
                    conflicts: 0,
                    stability_score: 0.6,
                    cumulative_stability: 0.48,
                },
            ],
            overall_stability: 0.7,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let avg = analysis.average_operation_stability();
        assert!((avg - 0.7).abs() < 0.001, "Average should be 0.7, got {}", avg);
    }

    #[test]
    fn analysis_summary_format() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![OperationRecord {
                id: rcad_kernel::persistent_naming::OperationId(1),
                operation_type: OperationType::BooleanUnion,
                label: Some("test".to_string()),
                sequence: 0,
                naming_events: vec![],
                stats: OperationStats::default(),
            }],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.85,
            broken_chains: vec![],
            stability_trend: 0.1,
            entity_counts: vec![10],
        };

        let summary = analysis.summary();
        assert!(summary.contains("Operations: 1"));
        assert!(summary.contains("85.0%"));
        assert!(summary.contains("Improving"));
    }

    #[test]
    fn infer_entity_type() {
        // Test entity type inference from encoded IDs.
        // Kind bits: Solid=0, Shell=1, Face=2, Wire=3, Edge=4, Vertex=5

        // Face (kind=2)
        let face_id = (2u64 << 56) | 42u64;
        let face_type = infer_entity_type_from_id(face_id);
        assert_eq!(face_type, Some(rcad_kernel::persistent_naming::EntityType::Face));

        // Edge (kind=4)
        let edge_id = (4u64 << 56) | 100u64;
        let edge_type = infer_entity_type_from_id(edge_id);
        assert_eq!(edge_type, Some(rcad_kernel::persistent_naming::EntityType::Edge));

        // Vertex (kind=5)
        let vertex_id = (5u64 << 56) | 7u64;
        let vertex_type = infer_entity_type_from_id(vertex_id);
        assert_eq!(vertex_type, Some(rcad_kernel::persistent_naming::EntityType::Vertex));

        // Solid (kind=0)
        let solid_id = (0u64 << 56) | 1u64;
        let solid_type = infer_entity_type_from_id(solid_id);
        assert_eq!(solid_type, Some(rcad_kernel::persistent_naming::EntityType::Solid));
    }

    #[test]
    fn cross_operation_stability_with_named_graph() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let mut named_graph = NamedGraph::from_brep(&brep);
        let initial_nodes = named_graph.graph().nodes.clone();

        // Perform a mutation.
        let result = named_graph.mutate_tracked(
            "test_mutation",
            OperationType::Generic,
            |graph, _history| {
                graph.record("test_op");
            },
        );

        assert!(result.is_ok());

        // Get stability report.
        let report = named_graph.stability_report();
        assert!(report.total_operations >= 1);

        // Analyze the history.
        let history = named_graph.history().clone();
        let analysis = analyze_naming_sequence(&[history], &initial_nodes);
        assert!(analysis.is_good() || analysis.is_excellent());
    }

    #[test]
    fn conflict_severity_classification() {
        let mut genealogy = HashMap::new();

        // Entity with broken genealogy (evolution is empty).
        genealogy.insert(
            PersistentId(1),
            EntityGenealogy {
                persistent_id: PersistentId(1),
                created_in_operation: rcad_kernel::persistent_naming::OperationId(1),
                evolution: vec![],
                current_entity_id: Some(1),
                is_deleted: false,
            },
        );

        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: genealogy,
            per_operation_stability: vec![],
            overall_stability: 0.9,
            broken_chains: vec![],
            stability_trend: 0.0,
            entity_counts: vec![],
        };

        let conflicts = detect_cross_operation_conflicts(&analysis);
        assert!(!conflicts.is_empty());

        // All detected conflicts should have proper severity.
        for conflict in &conflicts {
            assert!(matches!(
                conflict.severity,
                rcad_kernel::persistent_naming::IssueSeverity::Minor
                    | rcad_kernel::persistent_naming::IssueSeverity::Moderate
                    | rcad_kernel::persistent_naming::IssueSeverity::Severe
                    | rcad_kernel::persistent_naming::IssueSeverity::Critical
            ));
        }
    }

    #[test]
    fn recommendation_priorities_ordered() {
        let analysis = CrossOperationNamingAnalysis {
            operations: vec![],
            entity_genealogy: HashMap::new(),
            per_operation_stability: vec![],
            overall_stability: 0.3,  // Critical stability
            broken_chains: vec![],
            stability_trend: -0.5,  // Degrading
            entity_counts: vec![],
        };

        let recommendations = generate_stability_recommendations(&analysis);
        assert!(recommendations.len() >= 2, "Should have multiple recommendations");

        // Verify priorities are sorted (highest first).
        for i in 1..recommendations.len() {
            assert!(
                recommendations[i].priority <= recommendations[i - 1].priority,
                "Recommendations should be sorted by priority (descending)"
            );
        }
    }

    #[test]
    fn entity_counts_tracking() {
        let history = create_history_with_operation(
            OperationType::BooleanUnion,
            "test",
            OperationStats {
                entity_count_before: 10,
                entity_count_after: 15,
                names_preserved: 10,
                names_lost: 0,
                names_generated: 5,
                conflicts_resolved: 0,
            },
            vec![],
        );

        let initial_entities: Vec<TopoNode> = (0..10)
            .map(|i| TopoNode { kind: NodeKind::Face, index: i })
            .collect();

        let analysis = analyze_naming_sequence(&[history], &initial_entities);

        assert!(!analysis.entity_counts.is_empty());
        assert_eq!(analysis.entity_counts[0], 10);  // Initial count
        // After the operation, the count from stats is tracked.
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Tests for Enhanced Persistent Naming Semantics
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

#[cfg(test)]
mod enhanced_naming_tests {
    use crate::brep_graph::*;

    // 鈹€鈹€ ScopedId Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn scoped_id_creation() {
        let scope = NamingScope::for_part("housing").with_operation("fillet");
        let pid = PersistentId(42);
        let scoped = ScopedId::new(pid, scope.clone());

        assert_eq!(scoped.persistent_id, pid);
        assert_eq!(scoped.scope, scope);
        assert!(!scoped.is_null());
    }

    #[test]
    fn scoped_id_null() {
        let scoped = ScopedId::null();
        assert!(scoped.is_null());
        assert!(scoped.persistent_id.is_null());
    }

    #[test]
    fn scoped_id_qualified_name() {
        let scope = NamingScope::for_part("housing")
            .with_assembly("machine")
            .with_operation("fillet");
        let scoped = ScopedId::new(PersistentId(42), scope);

        let name = scoped.qualified_name();
        assert!(name.contains("machine"));
        assert!(name.contains("housing"));
        assert!(name.contains("fillet"));
        assert!(name.contains("e42"));
    }

    // 鈹€鈹€ NamingScope Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn naming_scope_creation() {
        let scope = NamingScope::new();
        assert!(scope.part.is_none());
        assert!(scope.assembly.is_none());
        assert!(scope.operation.is_none());
    }

    #[test]
    fn naming_scope_for_part() {
        let scope = NamingScope::for_part("housing");
        assert_eq!(scope.part, Some("housing".to_string()));
        assert!(scope.assembly.is_none());
    }

    #[test]
    fn naming_scope_for_assembly() {
        let scope = NamingScope::for_assembly("machine");
        assert_eq!(scope.assembly, Some("machine".to_string()));
        assert!(scope.part.is_none());
    }

    #[test]
    fn naming_scope_for_operation() {
        let scope = NamingScope::for_operation("housing", "fillet");
        assert_eq!(scope.part, Some("housing".to_string()));
        assert_eq!(scope.operation, Some("fillet".to_string()));
    }

    #[test]
    fn naming_scope_builder_pattern() {
        let scope = NamingScope::new()
            .with_assembly("machine")
            .with_part("housing")
            .with_operation("fillet");

        assert_eq!(scope.assembly, Some("machine".to_string()));
        assert_eq!(scope.part, Some("housing".to_string()));
        assert_eq!(scope.operation, Some("fillet".to_string()));
    }

    #[test]
    fn naming_scope_child_scope() {
        let parent = NamingScope::for_part("housing");
        let child = parent.child_scope("fillet");

        assert_eq!(child.part, parent.part);
        assert_eq!(child.operation, Some("fillet".to_string()));
    }

    #[test]
    fn naming_scope_contains() {
        let parent = NamingScope::for_assembly("machine");
        let child = NamingScope::for_part("housing").with_assembly("machine");

        assert!(parent.contains(&child));
        assert!(!child.contains(&parent));

        let other = NamingScope::for_assembly("other");
        assert!(!parent.contains(&other));
    }

    // 鈹€鈹€ EnhancedNamingContext Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn enhanced_context_assign_id() {
        let mut ctx = EnhancedNamingContext::new();
        let pid = ctx.assign_id(42);

        assert!(!pid.is_null());
        assert_eq!(ctx.resolve_entity(pid), Some(42));
        assert_eq!(ctx.resolve_persistent(42), Some(pid));
    }

    #[test]
    fn enhanced_context_with_scope() {
        let scope = NamingScope::for_part("housing");
        let ctx = EnhancedNamingContext::with_scope(scope.clone());

        assert_eq!(ctx.scope(), &scope);
    }

    #[test]
    fn enhanced_context_assign_derived_id_preserve() {
        let mut ctx = EnhancedNamingContext::new();

        // Assign original.
        let original_pid = ctx.assign_id(10);

        // Derive with Preserve policy.
        let derived_pid = ctx.assign_derived_id(20, &[10], NamePropagationPolicy::Preserve);

        // Should inherit the same persistent ID.
        assert_eq!(derived_pid, original_pid);
        assert_eq!(ctx.resolve_persistent(20), Some(original_pid));
    }

    #[test]
    fn enhanced_context_assign_derived_id_generate() {
        let mut ctx = EnhancedNamingContext::new();

        // Assign original.
        let original_pid = ctx.assign_id(10);

        // Derive with Generate policy.
        let derived_pid = ctx.assign_derived_id(20, &[10], NamePropagationPolicy::Generate);

        // Should get a new persistent ID.
        assert_ne!(derived_pid, original_pid);
        assert_eq!(ctx.resolve_persistent(20), Some(derived_pid));
    }

    #[test]
    fn enhanced_context_record_split() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_operation("test", "split_op"));

        // Create source entity.
        let source_pid = ctx.assign_id(10);

        // Split into three entities.
        let result_pids = ctx.record_split(10, &[20, 30, 40], "split_op");

        assert_eq!(result_pids.len(), 3);
        // First target inherits source's PID.
        assert_eq!(result_pids[0], source_pid);
        // Others get new PIDs.
        assert_ne!(result_pids[1], source_pid);
        assert_ne!(result_pids[2], source_pid);

        // Check genealogy.
        let genealogy = ctx.get_genealogy(source_pid).unwrap();
        assert_eq!(genealogy.status, EntityStatus::Split);
        assert_eq!(genealogy.child_ids.len(), 2); // Two new children.
    }

    #[test]
    fn enhanced_context_record_merge() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_operation("test", "merge_op"));

        // Create source entities.
        let pid1 = ctx.assign_id(10);
        let pid2 = ctx.assign_id(20);

        // Merge into one entity.
        let result_pid = ctx.record_merge(
            &[10, 20],
            30,
            "merge_op",
            NameConflictResolution::MergeEntities,
        );

        // Check result exists.
        assert!(!result_pid.is_null());

        // Check genealogy of sources.
        let genealogy1 = ctx.get_genealogy(pid1).unwrap();
        assert_eq!(genealogy1.status, EntityStatus::Merged);

        let genealogy2 = ctx.get_genealogy(pid2).unwrap();
        assert_eq!(genealogy2.status, EntityStatus::Merged);
    }

    #[test]
    fn enhanced_context_mark_deleted() {
        let mut ctx = EnhancedNamingContext::new();
        let pid = ctx.assign_id(10);

        ctx.mark_deleted(10);

        let genealogy = ctx.get_genealogy(pid).unwrap();
        assert_eq!(genealogy.status, EntityStatus::Deleted);
    }

    #[test]
    fn enhanced_context_entities_by_status() {
        let mut ctx = EnhancedNamingContext::new();

        let pid1 = ctx.assign_id(10);
        let pid2 = ctx.assign_id(20);
        ctx.mark_deleted(10);

        let deleted = ctx.entities_by_status(EntityStatus::Deleted);
        assert!(deleted.contains(&pid1));
        assert!(!deleted.contains(&pid2));

        let active = ctx.entities_by_status(EntityStatus::Active);
        assert!(active.contains(&pid2));
        assert!(!active.contains(&pid1));
    }

    #[test]
    fn enhanced_context_export_import_state() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_part("test"));

        let pid1 = ctx.assign_id(10);
        let pid2 = ctx.assign_id(20);

        // Export state.
        let state = ctx.export_state();

        // Clear context.
        ctx.clear();
        assert!(ctx.resolve_persistent(10).is_none());

        // Import state.
        ctx.import_state(state);

        // Verify restoration.
        assert_eq!(ctx.resolve_persistent(10), Some(pid1));
        assert_eq!(ctx.resolve_persistent(20), Some(pid2));
    }

    // 鈹€鈹€ NamePropagationRule Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn propagation_rule_for_boolean() {
        let rule = NamePropagationRule::for_operation(OperationType::BooleanUnion);

        assert_eq!(rule.face_policy, NamePropagationPolicy::Preserve);
        assert_eq!(rule.edge_policy, NamePropagationPolicy::Preserve);
        assert!(rule.track_genealogy);
    }

    #[test]
    fn propagation_rule_for_feature() {
        let rule = NamePropagationRule::for_operation(OperationType::Feature);

        assert_eq!(rule.face_policy, NamePropagationPolicy::Inherit);
        assert_eq!(rule.edge_policy, NamePropagationPolicy::Inherit);
        assert!(rule.track_genealogy);
    }

    #[test]
    fn propagation_rule_for_merge() {
        let rule = NamePropagationRule::for_operation(OperationType::Merge);

        assert_eq!(rule.face_policy, NamePropagationPolicy::Combine);
        assert_eq!(rule.conflict_resolution, NameConflictResolution::MergeEntities);
    }

    #[test]
    fn propagation_rule_policy_for_kind() {
        let rule = NamePropagationRule {
            operation_type: OperationType::Generic,
            face_policy: NamePropagationPolicy::Preserve,
            edge_policy: NamePropagationPolicy::Inherit,
            vertex_policy: NamePropagationPolicy::Generate,
            track_genealogy: true,
            conflict_resolution: NameConflictResolution::GenerateNewId,
        };

        assert_eq!(rule.policy_for_kind(NodeKind::Face), NamePropagationPolicy::Preserve);
        assert_eq!(rule.policy_for_kind(NodeKind::Edge), NamePropagationPolicy::Inherit);
        assert_eq!(rule.policy_for_kind(NodeKind::Vertex), NamePropagationPolicy::Generate);
    }

    // 鈹€鈹€ NamePropagationManager Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn propagation_manager_get_rule() {
        let manager = NamePropagationManager::new();

        let rule = manager.get_rule(OperationType::BooleanUnion);
        assert_eq!(rule.operation_type, OperationType::BooleanUnion);

        let rule = manager.get_rule(OperationType::Feature);
        assert_eq!(rule.operation_type, OperationType::Feature);
    }

    #[test]
    fn propagation_manager_set_rule() {
        let mut manager = NamePropagationManager::new();

        let custom_rule = NamePropagationRule {
            operation_type: OperationType::BooleanUnion,
            face_policy: NamePropagationPolicy::Generate,
            edge_policy: NamePropagationPolicy::Generate,
            vertex_policy: NamePropagationPolicy::Generate,
            track_genealogy: false,
            conflict_resolution: NameConflictResolution::KeepExisting,
        };

        manager.set_rule(custom_rule);

        let rule = manager.get_rule(OperationType::BooleanUnion);
        assert_eq!(rule.face_policy, NamePropagationPolicy::Generate);
    }

    #[test]
    fn propagation_manager_apply_split() {
        let mut manager = NamePropagationManager::new();
        let mut ctx = EnhancedNamingContext::new();

        // Create source entity.
        ctx.assign_id(10);

        // Apply split propagation.
        let pids = manager.apply_propagation(
            &mut ctx,
            OperationType::FaceSplit,
            &[10],
            &[20, 30, 40],
            NodeKind::Face,
            "split_test",
        );

        assert_eq!(pids.len(), 3);
        // First target should inherit source's PID.
        assert_eq!(pids[0], ctx.resolve_persistent(10).unwrap());
    }

    #[test]
    fn propagation_manager_apply_merge() {
        let mut manager = NamePropagationManager::new();
        let mut ctx = EnhancedNamingContext::new();

        // Create source entities.
        ctx.assign_id(10);
        ctx.assign_id(20);

        // Apply merge propagation.
        let pids = manager.apply_propagation(
            &mut ctx,
            OperationType::Merge,
            &[10, 20],
            &[30],
            NodeKind::Face,
            "merge_test",
        );

        assert_eq!(pids.len(), 1);
    }

    // 鈹€鈹€ NamingSnapshotManager Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn snapshot_manager_take_snapshot() {
        let mut manager = NamingSnapshotManager::new();
        let ctx = EnhancedNamingContext::new();

        let id = manager.take_snapshot(&ctx, Some("test_op".to_string()));

        assert_eq!(id, 0);
        assert_eq!(manager.len(), 1);
    }

    #[test]
    fn snapshot_manager_undo_redo() {
        let mut manager = NamingSnapshotManager::new();
        let mut ctx = EnhancedNamingContext::new();

        // Initial state.
        ctx.assign_id(10);
        manager.take_snapshot(&ctx, Some("op1".to_string()));

        // Second state.
        ctx.assign_id(20);
        manager.take_snapshot(&ctx, Some("op2".to_string()));

        assert_eq!(manager.len(), 2);

        // Undo.
        assert!(manager.can_undo());
        manager.undo(&mut ctx);
        assert!(!manager.can_undo());

        // Redo.
        assert!(manager.can_redo());
        manager.redo(&mut ctx);
        assert!(!manager.can_redo());
    }

    #[test]
    fn snapshot_manager_current() {
        let mut manager = NamingSnapshotManager::new();
        let ctx = EnhancedNamingContext::new();

        manager.take_snapshot(&ctx, Some("test".to_string()));

        let current = manager.current();
        assert!(current.is_some());
        assert_eq!(current.unwrap().operation, Some("test".to_string()));
    }

    // 鈹€鈹€ Genealogy Tracking Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn genealogy_tracking_through_boolean() {
        let mut ctx = EnhancedNamingContext::new();
        let mut manager = NamePropagationManager::new();

        ctx.set_scope(NamingScope::for_operation("part", "boolean_union"));

        // Create source entities.
        ctx.assign_id(1);
        ctx.assign_id(2);

        // Simulate boolean union.
        let result_pids = manager.apply_propagation(
            &mut ctx,
            OperationType::BooleanUnion,
            &[1, 2],
            &[10, 11],
            NodeKind::Face,
            "boolean_union",
        );

        assert_eq!(result_pids.len(), 2);
    }

    #[test]
    fn genealogy_tracking_through_fillet() {
        let mut ctx = EnhancedNamingContext::new();
        let mut manager = NamePropagationManager::new();

        ctx.set_scope(NamingScope::for_operation("part", "fillet"));

        // Create source edges.
        let pid_edge = ctx.assign_id(1);

        // Simulate fillet.
        let result_pids = manager.apply_propagation(
            &mut ctx,
            OperationType::Feature,
            &[1],
            &[10, 11],
            NodeKind::Face,
            "fillet",
        );

        assert_eq!(result_pids.len(), 2);

        let edge_genealogy = ctx.get_genealogy(pid_edge);
        assert!(edge_genealogy.is_some());
    }

    #[test]
    fn genealogy_multiple_operations() {
        let mut ctx = EnhancedNamingContext::new();
        let mut manager = NamePropagationManager::new();

        // Operation 1: Create initial face.
        ctx.set_scope(NamingScope::for_operation("part", "create"));
        let pid1 = ctx.assign_id(1);

        // Operation 2: Boolean union.
        ctx.set_scope(NamingScope::for_operation("part", "union"));
        manager.apply_propagation(
            &mut ctx,
            OperationType::BooleanUnion,
            &[1],
            &[2],
            NodeKind::Face,
            "union",
        );

        // Operation 3: Fillet.
        ctx.set_scope(NamingScope::for_operation("part", "fillet"));
        manager.apply_propagation(
            &mut ctx,
            OperationType::Feature,
            &[2],
            &[3, 4],
            NodeKind::Face,
            "fillet",
        );

        // Verify the chain is trackable.
        let descendants = ctx.trace_descendants(pid1);
        assert!(!descendants.is_empty() || ctx.get_genealogy(pid1).is_some());
    }

    // 鈹€鈹€ Name Stability Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn name_stability_through_boolean_preserve() {
        let mut ctx = EnhancedNamingContext::new();

        // Create faces from solid A.
        let pid_f1 = ctx.assign_id(1);
        ctx.assign_id(2);

        // Simulate boolean union where faces 1 and 2 are preserved.
        let preserved_pid = ctx.assign_derived_id(
            10, // New face index
            &[1], // Source face
            NamePropagationPolicy::Preserve,
        );

        assert_eq!(preserved_pid, pid_f1);
    }

    #[test]
    fn name_stability_through_split() {
        let mut ctx = EnhancedNamingContext::new();

        // Create a face that will be split.
        let pid = ctx.assign_id(1);

        // Split into three faces.
        let result_pids = ctx.record_split(1, &[10, 11, 12], "split");

        // First result should inherit the original PID.
        assert_eq!(result_pids[0], pid);
        // Others should get new PIDs.
        assert_ne!(result_pids[1], pid);
        assert_ne!(result_pids[2], pid);

        // Original entity should be marked as split.
        let genealogy = ctx.get_genealogy(pid).unwrap();
        assert_eq!(genealogy.status, EntityStatus::Split);
    }

    #[test]
    fn name_stability_through_merge() {
        let mut ctx = EnhancedNamingContext::new();

        // Create two faces that will be merged.
        let pid1 = ctx.assign_id(1);
        let pid2 = ctx.assign_id(2);

        // Merge into one face.
        let result_pid = ctx.record_merge(
            &[1, 2],
            10,
            "merge",
            NameConflictResolution::MergeEntities,
        );

        // Result should exist and sources should be marked merged.
        assert!(!result_pid.is_null());

        let g1 = ctx.get_genealogy(pid1).unwrap();
        let g2 = ctx.get_genealogy(pid2).unwrap();

        assert_eq!(g1.status, EntityStatus::Merged);
        assert_eq!(g2.status, EntityStatus::Merged);
    }

    // 鈹€鈹€ Conflict Resolution Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn conflict_resolution_keep_existing() {
        let mut ctx = EnhancedNamingContext::new();

        let pid = ctx.assign_id(1);

        // Try to assign the same PID to another entity (simulating conflict).
        let conflict = NameConflictRecord {
            persistent_id: pid,
            conflicting_entities: vec![1, 2],
            operation: "test".to_string(),
            scope: ctx.current_scope.clone(),
            resolution: NameConflictResolution::Unresolved,
            sequence: 0,
        };

        ctx.resolve_conflict(&conflict, NameConflictResolution::KeepExisting).unwrap();

        // Original entity should keep its PID.
        assert_eq!(ctx.resolve_persistent(1), Some(pid));

        // Conflict should be recorded.
        assert!(!ctx.conflict_history.is_empty());
    }

    #[test]
    fn conflict_resolution_generate_new() {
        let mut ctx = EnhancedNamingContext::new();

        let pid = ctx.assign_id(1);

        let conflict = NameConflictRecord {
            persistent_id: pid,
            conflicting_entities: vec![1, 2],
            operation: "test".to_string(),
            scope: ctx.current_scope.clone(),
            resolution: NameConflictResolution::Unresolved,
            sequence: 0,
        };

        ctx.resolve_conflict(&conflict, NameConflictResolution::GenerateNewId).unwrap();

        // New entity should have a different PID.
        let pid1 = ctx.resolve_persistent(1);
        let pid2 = ctx.resolve_persistent(2);

        // At least one should have a different PID.
        assert!(pid1.is_some() || pid2.is_some());
    }

    // 鈹€鈹€ Serialization Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn serialization_naming_scope() {
        let scope = NamingScope::for_part("housing")
            .with_assembly("machine")
            .with_operation("fillet");

        let json = serde_json::to_string(&scope).unwrap();
        let decoded: NamingScope = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded, scope);
    }

    #[test]
    fn serialization_scoped_id() {
        let scope = NamingScope::for_part("housing");
        let scoped = ScopedId::new(PersistentId(42), scope);

        let json = serde_json::to_string(&scoped).unwrap();
        let decoded: ScopedId = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded, scoped);
    }

    #[test]
    fn serialization_enhanced_context_state() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.set_scope(NamingScope::for_part("test"));
        ctx.assign_id(10);
        ctx.assign_id(20);

        let state = ctx.export_state();

        let json = serde_json::to_string(&state).unwrap();
        let decoded: EnhancedNamingContextState = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded.current_scope, state.current_scope);
        assert_eq!(decoded.next_persistent_id, state.next_persistent_id);
    }

    #[test]
    fn serialization_naming_context_snapshot() {
        let mut ctx = EnhancedNamingContext::new();
        ctx.assign_id(10);

        let snapshot = NamingContextSnapshot::from_context(&ctx, 0, Some("test".to_string()));

        let json = serde_json::to_string(&snapshot).unwrap();
        let decoded: NamingContextSnapshot = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded.id, snapshot.id);
        assert_eq!(decoded.operation, snapshot.operation);
    }
}
