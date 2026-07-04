#[cfg(test)]
mod tests {
    use crate::history::*;
    use rcad_kernel::{PrimitiveSolid, TopoEntityRef};

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    // 鈹€鈹€ HistoryTracker Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn tracker_starts_empty() {
        let tracker = HistoryTracker::new();
        assert!(!tracker.has_modified());
        assert!(!tracker.has_generated());
        assert!(!tracker.has_deleted());
    }

    #[test]
    fn tracker_records_face_modification() {
        let mut tracker = HistoryTracker::new();
        tracker.record_face_modified(0, 1);

        assert!(tracker.has_modified());
        assert!(tracker.is_face_modified(0));
        assert_eq!(tracker.modified_faces(0), vec![1]);

        let source = tracker.get_source(EntityType::Face, 1);
        assert_eq!(source, Some((InputSource::A, 0)));
    }

    #[test]
    fn tracker_records_face_split() {
        let mut tracker = HistoryTracker::new();
        tracker.record_face_modified_multi(0, vec![1, 2, 3], InputSource::A);

        assert!(tracker.has_modified());
        let record = tracker.face_modification_record(0).unwrap();
        assert_eq!(record.modification_type, ModificationType::Split);
        assert_eq!(record.result_indices, vec![1, 2, 3]);
    }

    #[test]
    fn tracker_records_generated_entities() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_generated(10, GenerationCause::Intersection);
        tracker.record_edge_generated(20, GenerationCause::NewBoundary);
        tracker.record_vertex_generated(30, GenerationCause::Intersection);

        assert!(tracker.has_generated());
        assert!(tracker.is_face_generated(10));
        assert!(tracker.is_edge_generated(20));
        assert!(tracker.is_vertex_generated(30));

        assert!(!tracker.is_face_generated(11));
    }

    #[test]
    fn tracker_records_deletions() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_deleted(0, DeletionReason::BooleanOperation);
        tracker.record_edge_deleted(1, DeletionReason::Overlap);
        tracker.record_vertex_deleted(2, DeletionReason::Custom("Test".to_string()));

        assert!(tracker.has_deleted());
        assert!(tracker.is_face_deleted(0));
        assert!(tracker.is_edge_deleted(1));
        assert!(tracker.is_vertex_deleted(2));

        let record = tracker.deletion_record(0, EntityType::Face).unwrap();
        assert_eq!(record.reason, DeletionReason::BooleanOperation);

        let vertex_record = tracker.deletion_record(2, EntityType::Vertex).unwrap();
        if let DeletionReason::Custom(s) = &vertex_record.reason {
            assert_eq!(s, "Test");
        } else {
            panic!("Expected Custom deletion reason");
        }
    }

    #[test]
    fn tracker_count_by_source() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_modified_multi(0, vec![1, 2], InputSource::A);
        tracker.record_face_modified_multi(1, vec![3], InputSource::B);
        tracker.record_face_generated(4, GenerationCause::Intersection);

        assert_eq!(tracker.count_by_source(EntityType::Face, InputSource::A), 2);
        assert_eq!(tracker.count_by_source(EntityType::Face, InputSource::B), 1);
        assert_eq!(tracker.count_by_source(EntityType::Face, InputSource::Generated), 1);
    }

    #[test]
    fn tracker_merge() {
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified(0, 1);

        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified(2, 3);
        tracker2.record_face_generated(4, GenerationCause::Intersection);

        tracker1.merge(&tracker2);

        assert!(tracker1.is_face_modified(0));
        assert!(tracker1.is_face_modified(2));
        assert!(tracker1.is_face_generated(4));
    }

    #[test]
    fn tracker_statistics() {
        let mut tracker = HistoryTracker::new();

        tracker.record_face_modified(0, 1);
        tracker.record_edge_modified(0, vec![1, 2], InputSource::A, ModificationType::Split);
        tracker.record_face_generated(10, GenerationCause::Intersection);
        tracker.record_face_deleted(20, DeletionReason::BooleanOperation);

        let stats = tracker.statistics();
        assert_eq!(stats.modified_faces, 1);
        assert_eq!(stats.modified_edges, 1);
        assert_eq!(stats.generated_faces, 1);
        assert_eq!(stats.deleted_faces, 1);
        assert_eq!(stats.total_modified(), 2);
        assert_eq!(stats.total_generated(), 1);
        assert_eq!(stats.total_deleted(), 1);
    }

    // 鈹€鈹€ BooleanHistory Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn boolean_history_new() {
        let history = BooleanHistory::new();
        assert!(history.is_empty());
        assert!(!history.has_modified());
        assert!(!history.has_generated());
        assert!(!history.has_deleted());
    }

    #[test]
    fn boolean_history_populate_tracker() {
        let mut history = BooleanHistory::new();
        history.face_origins = vec![
            FaceOrigin::FromA(0),
            FaceOrigin::FromB(1),
            FaceOrigin::Generated,
        ];
        history.edge_origins = vec![
            EdgeOrigin::FromA(0),
            EdgeOrigin::SplitFromA(0),
            EdgeOrigin::Generated,
        ];
        history.vertex_origins = vec![
            VertexOrigin::FromA(0),
            VertexOrigin::Intersection,
        ];

        history.populate_tracker();

        assert!(history.has_modified());
        assert!(history.has_generated());

        // Check face queries.
        assert_eq!(history.modified_faces(0, true), vec![0]);
        assert_eq!(history.modified_faces(1, false), vec![1]);
        assert_eq!(history.generated_faces(), vec![2]);

        // Check edge queries.
        let modified_edges = history.modified_edges(0, true);
        assert_eq!(modified_edges.len(), 2);
        assert!(modified_edges.contains(&0));

        // Check vertex queries.
        assert_eq!(history.generated_vertices(), vec![1]);
    }

    #[test]
    fn boolean_history_get_source() {
        let mut history = BooleanHistory::new();
        history.face_origins = vec![
            FaceOrigin::FromA(0),
            FaceOrigin::FromB(5),
            FaceOrigin::Generated,
        ];

        assert_eq!(history.get_face_source(0), Some((0, true)));
        assert_eq!(history.get_face_source(1), Some((5, false)));
        assert_eq!(history.get_face_source(2), None);
    }

    #[test]
    fn boolean_history_deletions() {
        let mut history = BooleanHistory::new();
        history.deleted_from_a = vec![0, 1];
        history.deleted_from_b = vec![2];
        history.deletion_reasons.insert(
            (EntityType::Face, 0),
            DeletionReason::BooleanOperation,
        );

        assert!(history.is_face_deleted(0, true));
        assert!(history.is_face_deleted(1, true));
        assert!(history.is_face_deleted(2, false));
        assert!(!history.is_face_deleted(3, true));
    }

    #[test]
    fn propagate_persistent_naming_maps_face_edge_vertex_and_solid_origins() {
        let result_brep = unit_box();
        let mut history = BooleanHistory {
            face_origins: vec![FaceOrigin::FromA(0), FaceOrigin::FromB(1)],
            co_face_origins: vec![],
            edge_origins: vec![EdgeOrigin::FromA(0), EdgeOrigin::SplitFromA(0), EdgeOrigin::FromB(1)],
            vertex_origins: vec![VertexOrigin::FromA(0), VertexOrigin::Intersection, VertexOrigin::FromB(1)],
            shell_origins: vec![],
            solid_origins: vec![SolidOrigin::FromA],
            tracker: HistoryTracker::new(),
            deleted_from_a: vec![],
            deleted_from_b: vec![],
            deletion_reasons: HashMap::new(),
            source_history: Vec::new(),
        };

        let mut names_a = PersistentNamingHooks::new();
        names_a.bind("face_a", TopoEntityRef::Face(0));
        names_a.bind("edge_a", TopoEntityRef::Edge(0));
        names_a.bind("vertex_a", TopoEntityRef::Vertex(0));
        names_a.bind("solid_a", TopoEntityRef::Solid(0));

        let mut names_b = PersistentNamingHooks::new();
        names_b.bind("face_b", TopoEntityRef::Face(1));
        names_b.bind("edge_b", TopoEntityRef::Edge(1));
        names_b.bind("vertex_b", TopoEntityRef::Vertex(1));
        names_b.bind("solid_b", TopoEntityRef::Solid(0));

        let (result_names, report) = history.propagate_persistent_naming(&result_brep, &names_a, &names_b);

        assert_eq!(result_names.resolve("face_a"), Some(TopoEntityRef::Face(0)));
        assert_eq!(result_names.resolve("face_b"), Some(TopoEntityRef::Face(1)));
        assert_eq!(result_names.resolve("edge_a"), Some(TopoEntityRef::Edge(0)));
        assert_eq!(result_names.resolve("edge_a@1"), Some(TopoEntityRef::Edge(1)));
        assert_eq!(result_names.resolve("edge_b"), Some(TopoEntityRef::Edge(2)));
        assert_eq!(result_names.resolve("vertex_a"), Some(TopoEntityRef::Vertex(0)));
        assert_eq!(result_names.resolve("vertex_b"), Some(TopoEntityRef::Vertex(2)));
        assert_eq!(result_names.resolve("solid_a"), Some(TopoEntityRef::Solid(0)));
        assert_eq!(result_names.resolve("solid_b"), None);

        assert!(report.dropped_from_a.is_empty());
        assert_eq!(report.dropped_from_b, vec!["solid_b".to_string()]);
        assert_eq!(report.duplicated_from_a, vec!["edge_a".to_string()]);
        assert!(report.duplicated_from_b.is_empty());
    }

    #[test]
    fn propagate_persistent_naming_disambiguates_cross_input_name_collisions() {
        let result_brep = unit_box();
        let mut history = BooleanHistory {
            face_origins: vec![FaceOrigin::FromA(0), FaceOrigin::FromB(0)],
            co_face_origins: vec![],
            edge_origins: vec![],
            vertex_origins: vec![],
            shell_origins: vec![],
            solid_origins: vec![],
            tracker: HistoryTracker::new(),
            deleted_from_a: vec![],
            deleted_from_b: vec![],
            deletion_reasons: HashMap::new(),
            source_history: Vec::new(),
        };        let mut names_a = PersistentNamingHooks::new();
        names_a.bind("shared_face", TopoEntityRef::Face(0));
        let mut names_b = PersistentNamingHooks::new();
        names_b.bind("shared_face", TopoEntityRef::Face(0));

        let (result_names, report) = history.propagate_persistent_naming(&result_brep, &names_a, &names_b);

        assert_eq!(result_names.resolve("shared_face"), Some(TopoEntityRef::Face(0)));
        assert_eq!(result_names.resolve("shared_face@1"), Some(TopoEntityRef::Face(1)));
        assert!(report.dropped_from_a.is_empty());
        assert!(report.dropped_from_b.is_empty());
    }

    // 鈹€鈹€ HistoryChain Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn chain_starts_empty() {
        let chain = HistoryChain::new();
        assert!(chain.is_empty());
        assert_eq!(chain.len(), 0);
    }

    #[test]
    fn chain_push_and_get() {
        let mut chain = HistoryChain::new();

        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified(0, 1);
        chain.push(tracker1, Some("op1".to_string()));

        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified(1, 2);
        chain.push(tracker2, Some("op2".to_string()));

        assert_eq!(chain.len(), 2);
        assert_eq!(chain.label(0), Some("op1"));
        assert_eq!(chain.label(1), Some("op2"));
        assert!(chain.get(0).unwrap().is_face_modified(0));
        assert!(chain.get(1).unwrap().is_face_modified(1));
    }

    #[test]
    fn chain_trace_ancestry() {
        let mut chain = HistoryChain::new();

        // First operation: face 0 -> face 1.
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified_multi(0, vec![1], InputSource::A);
        chain.push(tracker1, None);

        // Second operation: face 1 -> face 2.
        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified_multi(1, vec![2], InputSource::A);
        chain.push(tracker2, None);

        // Trace ancestry of face 2.
        let ancestry = chain.trace_ancestry(EntityType::Face, 2);
        assert_eq!(ancestry.len(), 2);

        // Should trace back through both operations.
        assert_eq!(ancestry[0].0, 1); // Second operation.
        assert_eq!(ancestry[0].1, 2); // Face 2.

        assert_eq!(ancestry[1].0, 0); // First operation.
        assert_eq!(ancestry[1].1, 1); // Face 1.
    }

    #[test]
    fn chain_trace_descendants() {
        let mut chain = HistoryChain::new();

        // First operation: face 0 -> faces 1, 2 (split).
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified_multi(0, vec![1, 2], InputSource::A);
        chain.push(tracker1, None);

        // Second operation: faces 1, 2 -> faces 3, 4, 5.
        let mut tracker2 = HistoryTracker::new();
        tracker2.record_face_modified_multi(1, vec![3], InputSource::A);
        tracker2.record_face_modified_multi(2, vec![4, 5], InputSource::A);
        chain.push(tracker2, None);

        // Trace descendants of face 0.
        let descendants = chain.trace_descendants(EntityType::Face, 0);
        assert_eq!(descendants.len(), 2);

        // First operation: 0 -> [1, 2].
        assert_eq!(descendants[0].0, 0);
        assert_eq!(descendants[0].1.len(), 2);

        // Second operation: 1, 2 -> [3, 4, 5].
        assert_eq!(descendants[1].0, 1);
        assert_eq!(descendants[1].1.len(), 3);
    }

    #[test]
    fn chain_is_deleted_any() {
        let mut chain = HistoryChain::new();

        // First operation: delete face 0.
        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_deleted(0, DeletionReason::BooleanOperation);
        chain.push(tracker1, None);

        // Second operation: no deletions.
        let tracker2 = HistoryTracker::new();
        chain.push(tracker2, None);

        assert!(chain.is_deleted_any(EntityType::Face, 0));
        assert!(!chain.is_deleted_any(EntityType::Face, 1));
    }

    #[test]
    fn chain_statistics() {
        let mut chain = HistoryChain::new();

        let mut tracker1 = HistoryTracker::new();
        tracker1.record_face_modified(0, 1);
        tracker1.record_face_generated(2, GenerationCause::Intersection);
        chain.push(tracker1, None);

        let mut tracker2 = HistoryTracker::new();
        tracker2.record_edge_modified(0, vec![1, 2], InputSource::A, ModificationType::Split);
        tracker2.record_face_deleted(3, DeletionReason::BooleanOperation);
        chain.push(tracker2, None);

        let stats = chain.statistics();
        assert_eq!(stats.operation_count, 2);
        assert_eq!(stats.total_modified_faces, 1);
        assert_eq!(stats.total_modified_edges, 1);
        assert_eq!(stats.total_generated_faces, 1);
        assert_eq!(stats.total_deleted_faces, 1);
    }

    // 鈹€鈹€ DeletionReason Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn deletion_reason_description() {
        assert_eq!(DeletionReason::BooleanOperation.description(), "Removed by boolean operation");
        assert_eq!(DeletionReason::OutsideResult.description(), "Outside result volume");
        assert_eq!(DeletionReason::Overlap.description(), "Overlapping geometry");
        assert_eq!(DeletionReason::Tolerance.description(), "Tolerance issues");
        assert_eq!(DeletionReason::Healing.description(), "Removed during healing");
        assert_eq!(
            DeletionReason::Custom("Custom reason".to_string()).description(),
            "Custom reason"
        );
    }

    // 鈹€鈹€ GenerationCause Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn generation_record_with_parents() {
        let mut tracker = HistoryTracker::new();
        tracker.record_edge_generated_with_parents(
            10,
            GenerationCause::Intersection,
            vec![0, 1], // Parent edges
        );

        let record = tracker.edge_generation_record(10).unwrap();
        assert_eq!(record.entity_index, 10);
        assert_eq!(record.cause, GenerationCause::Intersection);
        assert_eq!(record.parent_indices, vec![0, 1]);
    }
}
