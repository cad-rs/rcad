use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_kernel::persistent_naming::{
    PersistentId, PersistentNamingEngine, NamingStabilityReport,
    OperationType, OperationStats, CrossOperationStabilityReport,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NodeKind {
    Solid,
    Shell,
    Face,
    Wire,
    Edge,
    Vertex,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TopoNode {
    pub kind: NodeKind,
    pub index: usize,
}

#[derive(Debug, Clone)]
pub struct TopoGraphHistoryEvent {
    pub action: String,
}

#[derive(Debug, Clone, Default)]
pub struct TopoGraphHistory {
    pub events: Vec<TopoGraphHistoryEvent>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TopoGraphValidationIssue {
    MissingAdjacency { node: TopoNode },
    NonSymmetricAdjacency { a: TopoNode, b: TopoNode },
    InvalidEdgeVertexRef { edge_index: usize, vertex_index: usize },
}

#[derive(Debug, Clone, Default)]
pub struct TopoGraph {
    pub nodes: Vec<TopoNode>,
    pub history: TopoGraphHistory,
    adjacency: HashMap<TopoNode, Vec<TopoNode>>,
    solid_shells: Vec<Vec<usize>>,
    shell_faces: Vec<Vec<usize>>,
    face_wires: Vec<Vec<usize>>,
    wire_edges: Vec<Vec<usize>>,
    edge_vertices: Vec<[usize; 2]>,
}

impl TopoGraph {
    pub fn from_brep(brep: &BRep) -> Self {
        let mut g = Self::default();
        g.record("from_brep");

        for vi in 0..brep.vertices.len() {
            g.add_node(TopoNode {
                kind: NodeKind::Vertex,
                index: vi,
            });
        }

        for (ei, e) in brep.edges.iter().enumerate() {
            let en = TopoNode {
                kind: NodeKind::Edge,
                index: ei,
            };
            g.add_node(en);
            g.edge_vertices.push([e.start, e.end]);
            g.connect(
                en,
                TopoNode {
                    kind: NodeKind::Vertex,
                    index: e.start,
                },
            );
            g.connect(
                en,
                TopoNode {
                    kind: NodeKind::Vertex,
                    index: e.end,
                },
            );
        }

        let mut shell_idx = 0usize;
        let mut face_idx = 0usize;
        let mut wire_idx = 0usize;

        for (si, solid) in brep.solids.iter().enumerate() {
            let sn = TopoNode {
                kind: NodeKind::Solid,
                index: si,
            };
            g.add_node(sn);
            g.solid_shells.push(Vec::new());

            for shell in &solid.shells {
                let shn = TopoNode {
                    kind: NodeKind::Shell,
                    index: shell_idx,
                };
                g.add_node(shn);
                g.connect(sn, shn);
                g.solid_shells[si].push(shell_idx);
                g.shell_faces.push(Vec::new());

                for face in &shell.faces {
                    let fnn = TopoNode {
                        kind: NodeKind::Face,
                        index: face_idx,
                    };
                    g.add_node(fnn);
                    g.connect(shn, fnn);
                    g.shell_faces[shell_idx].push(face_idx);
                    g.face_wires.push(Vec::new());

                    let wires = std::iter::once(&face.outer_wire).chain(face.inner_wires.iter());
                    for wire in wires {
                        let wn = TopoNode {
                            kind: NodeKind::Wire,
                            index: wire_idx,
                        };
                        g.add_node(wn);
                        g.connect(fnn, wn);
                        g.face_wires[face_idx].push(wire_idx);
                        g.wire_edges.push(Vec::new());

                        for we in &wire.edges {
                            let en = TopoNode {
                                kind: NodeKind::Edge,
                                index: we.idx,
                            };
                            if we.idx < brep.edges.len() {
                                g.connect(wn, en);
                                if !g.wire_edges[wire_idx].contains(&we.idx) {
                                    g.wire_edges[wire_idx].push(we.idx);
                                }
                            }
                        }

                        wire_idx += 1;
                    }

                    face_idx += 1;
                }

                shell_idx += 1;
            }
        }

        g
    }

    pub fn record(&mut self, action: impl Into<String>) {
        self.history.events.push(TopoGraphHistoryEvent {
            action: action.into(),
        });
    }

    pub fn neighbors(&self, node: TopoNode) -> Vec<TopoNode> {
        self.adjacency.get(&node).cloned().unwrap_or_default()
    }

    pub fn faces_of_shell(&self, shell: TopoNode) -> Vec<TopoNode> {
        if shell.kind != NodeKind::Shell {
            return Vec::new();
        }
        self.shell_faces
            .get(shell.index)
            .map(|v| {
                v.iter()
                    .map(|&i| TopoNode {
                        kind: NodeKind::Face,
                        index: i,
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    pub fn edges_of_face(&self, face: TopoNode) -> Vec<TopoNode> {
        if face.kind != NodeKind::Face {
            return Vec::new();
        }
        let mut out: Vec<usize> = Vec::new();
        if let Some(wires) = self.face_wires.get(face.index) {
            for &wi in wires {
                if let Some(edges) = self.wire_edges.get(wi) {
                    for &ei in edges {
                        if !out.contains(&ei) {
                            out.push(ei);
                        }
                    }
                }
            }
        }
        out.into_iter()
            .map(|i| TopoNode {
                kind: NodeKind::Edge,
                index: i,
            })
            .collect()
    }

    pub fn vertices_of_edge(&self, edge: TopoNode) -> Vec<TopoNode> {
        if edge.kind != NodeKind::Edge {
            return Vec::new();
        }
        self.edge_vertices
            .get(edge.index)
            .map(|v| {
                vec![
                    TopoNode {
                        kind: NodeKind::Vertex,
                        index: v[0],
                    },
                    TopoNode {
                        kind: NodeKind::Vertex,
                        index: v[1],
                    },
                ]
            })
            .unwrap_or_default()
    }

    pub fn validate(&self) -> Vec<TopoGraphValidationIssue> {
        let mut issues = Vec::new();
        for node in &self.nodes {
            let Some(neigh) = self.adjacency.get(node) else {
                issues.push(TopoGraphValidationIssue::MissingAdjacency { node: *node });
                continue;
            };
            for n in neigh {
                if let Some(back) = self.adjacency.get(n) {
                    if !back.contains(node) {
                        issues.push(TopoGraphValidationIssue::NonSymmetricAdjacency {
                            a: *node,
                            b: *n,
                        });
                    }
                } else {
                    issues.push(TopoGraphValidationIssue::MissingAdjacency { node: *n });
                }
            }
        }

        for (ei, vv) in self.edge_vertices.iter().enumerate() {
            for &vi in vv {
                if !self.nodes.contains(&TopoNode {
                    kind: NodeKind::Vertex,
                    index: vi,
                }) {
                    issues.push(TopoGraphValidationIssue::InvalidEdgeVertexRef {
                        edge_index: ei,
                        vertex_index: vi,
                    });
                }
            }
        }

        issues
    }

    /// Compact graph storage by dropping orphan adjacency entries and
    /// deduplicating neighbor lists.
    pub fn compact(&mut self) {
        let mut node_set = std::collections::HashSet::new();
        for n in &self.nodes {
            node_set.insert(*n);
        }

        self.adjacency.retain(|node, _| node_set.contains(node));
        for neigh in self.adjacency.values_mut() {
            neigh.retain(|n| node_set.contains(n));
            neigh.sort_by_key(|n| {
                let kind_rank = match n.kind {
                    NodeKind::Solid => 0usize,
                    NodeKind::Shell => 1,
                    NodeKind::Face => 2,
                    NodeKind::Wire => 3,
                    NodeKind::Edge => 4,
                    NodeKind::Vertex => 5,
                };
                (kind_rank, n.index)
            });
            neigh.dedup();
        }
        self.record("compact");
    }

    /// Apply a mutation and run graph validation afterward.
    ///
    /// This is a lightweight baseline for mutation-guard workflows: callers
    /// can route all topology edits through this helper and reject invalid
    /// states before continuing downstream processing.
    pub fn mutate_checked<F>(
        &mut self,
        action: impl Into<String>,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph),
    {
        let action = action.into();
        mutator(self);
        let issues = self.validate();
        if issues.is_empty() {
            self.record(format!("mutate:{action}"));
            Ok(())
        } else {
            self.record(format!("mutate_invalid:{action}"));
            Err(issues)
        }
    }

    /// Apply a mutation with rollback-on-failure semantics.
    ///
    /// If validation fails after the mutation, graph state is restored to the
    /// pre-mutation snapshot and validation issues are returned.
    pub fn mutate_guarded<F>(
        &mut self,
        action: impl Into<String>,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph),
    {
        let action = action.into();
        let before = self.clone();
        mutator(self);
        let issues = self.validate();
        if issues.is_empty() {
            self.record(format!("mutate_guarded:{action}"));
            Ok(())
        } else {
            *self = before;
            self.record(format!("mutate_guarded_rollback:{action}"));
            Err(issues)
        }
    }

    fn add_node(&mut self, node: TopoNode) {
        if !self.nodes.contains(&node) {
            self.nodes.push(node);
        }
        self.adjacency.entry(node).or_default();
    }

    fn connect(&mut self, a: TopoNode, b: TopoNode) {
        self.add_node(a);
        self.add_node(b);
        let va = self.adjacency.entry(a).or_default();
        if !va.contains(&b) {
            va.push(b);
        }
        let vb = self.adjacency.entry(b).or_default();
        if !vb.contains(&a) {
            vb.push(a);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRepGraphHistory: Persistent Naming Integration
// ─────────────────────────────────────────────────────────────────────────────

/// Enhanced history with persistent naming integration for cross-operation stability.
///
/// This struct bridges the TopoGraph mutation history with the PersistentNamingEngine,
/// enabling:
/// - Automatic name propagation during topology mutations
/// - Cross-operation stability analysis
/// - Entity genealogy tracking
/// - Undo/redo support with naming reconstruction
#[derive(Debug, Clone)]
pub struct BRepGraphHistory {
    /// The underlying naming engine.
    naming_engine: PersistentNamingEngine,
    /// Snapshots for undo support.
    snapshots: Vec<TopoGraphSnapshot>,
    /// Current snapshot index (for undo/redo).
    current_snapshot: usize,
}

/// A snapshot of the graph state with naming context.
#[derive(Debug, Clone)]
struct TopoGraphSnapshot {
    /// The action that created this snapshot.
    action: String,
    /// Node count at this snapshot.
    node_count: usize,
    /// Entity ID to persistent ID mappings at this snapshot.
    naming: HashMap<TopoNode, PersistentId>,
}

impl Default for BRepGraphHistory {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepGraphHistory {
    /// Create a new history with default naming engine.
    pub fn new() -> Self {
        Self {
            naming_engine: PersistentNamingEngine::default(),
            snapshots: Vec::new(),
            current_snapshot: 0,
        }
    }

    /// Create a history with a specific naming rule.
    pub fn with_naming_rule(rule: rcad_kernel::persistent_naming::NamingRule) -> Self {
        Self {
            naming_engine: PersistentNamingEngine::new(rule),
            snapshots: Vec::new(),
            current_snapshot: 0,
        }
    }

    /// Get a reference to the naming engine.
    pub fn naming_engine(&self) -> &PersistentNamingEngine {
        &self.naming_engine
    }

    /// Get mutable access to the naming engine.
    pub fn naming_engine_mut(&mut self) -> &mut PersistentNamingEngine {
        &mut self.naming_engine
    }

    /// Begin an operation for history tracking.
    pub fn begin_operation(&mut self, operation_type: OperationType, label: Option<String>) {
        self.naming_engine.begin_operation(operation_type, label);
    }

    /// Finalize the current operation.
    pub fn finalize_operation(&mut self, stats: OperationStats) {
        self.naming_engine.finalize_operation(stats);
    }

    /// Record a graph mutation with naming propagation.
    ///
    /// This creates a snapshot and propagates names for surviving entities.
    pub fn record_mutation(
        &mut self,
        graph: &TopoGraph,
        action: &str,
        entity_map: &[(TopoNode, Option<TopoNode>)],
    ) {
        // Create snapshot.
        let mut naming = HashMap::new();
        for node in &graph.nodes {
            if let Some(pid) = self.naming_engine.resolve_persistent(node_to_entity_id(*node)) {
                naming.insert(*node, pid);
            }
        }

        self.snapshots.truncate(self.current_snapshot + 1);
        self.snapshots.push(TopoGraphSnapshot {
            action: action.to_string(),
            node_count: graph.nodes.len(),
            naming,
        });
        self.current_snapshot = self.snapshots.len() - 1;

        // Propagate names for surviving entities.
        let entity_id_map: Vec<(u64, Option<u64>)> = entity_map
            .iter()
            .map(|(old, new)| {
                (
                    node_to_entity_id(*old),
                    new.map(node_to_entity_id),
                )
            })
            .collect();

        self.naming_engine.propagate_names(
            &entity_id_map,
            rcad_kernel::persistent_naming::NamePropagationPolicy::Preserve,
        );
    }

    /// Assign a persistent ID to a topology node.
    pub fn assign_persistent_id(&mut self, node: TopoNode) -> PersistentId {
        self.naming_engine.assign_persistent_id(node_to_entity_id(node))
    }

    /// Resolve a topology node to its persistent ID.
    pub fn resolve_persistent(&self, node: TopoNode) -> Option<PersistentId> {
        self.naming_engine.resolve_persistent(node_to_entity_id(node))
    }

    /// Resolve a persistent ID back to a topology node.
    pub fn resolve_node(&self, pid: PersistentId, kind: NodeKind) -> Option<TopoNode> {
        let entity_id = self.naming_engine.resolve_entity(pid)?;
        Some(TopoNode {
            kind,
            index: entity_id as usize,
        })
    }

    /// Generate a stability report for the current state.
    pub fn stability_report(&self) -> CrossOperationStabilityReport {
        self.naming_engine.cross_operation_stability_report()
    }

    /// Generate a naming stability report comparing before and after states.
    pub fn naming_stability_report(
        &self,
        _before_nodes: &[TopoNode],
        after_nodes: &[TopoNode],
    ) -> NamingStabilityReport {
        let before_context = self.naming_engine.context().clone();
        let after_ids: Vec<u64> = after_nodes.iter().map(|n| node_to_entity_id(*n)).collect();
        self.naming_engine.stability_report(&before_context, &after_ids)
    }

    /// Track an edge split event.
    pub fn track_edge_split(
        &mut self,
        old_edge_idx: usize,
        new_edge_indices: &[usize],
    ) -> Vec<PersistentId> {
        self.naming_engine.propagate_split(
            old_edge_idx as u64,
            &new_edge_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
        )
    }

    /// Track a face split event.
    pub fn track_face_split(
        &mut self,
        old_face_idx: usize,
        new_face_indices: &[usize],
    ) -> Vec<PersistentId> {
        self.naming_engine.propagate_split(
            old_face_idx as u64,
            &new_face_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
        )
    }

    /// Track a vertex merge event.
    pub fn track_vertex_merge(
        &mut self,
        old_vertex_indices: &[usize],
        new_vertex_idx: usize,
    ) -> PersistentId {
        self.naming_engine.propagate_merge(
            &old_vertex_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
            new_vertex_idx as u64,
            rcad_kernel::persistent_naming::NamePropagationPolicy::Preserve,
        )
    }

    /// Track a face merge event.
    pub fn track_face_merge(
        &mut self,
        old_face_indices: &[usize],
        new_face_idx: usize,
    ) -> PersistentId {
        self.naming_engine.propagate_merge(
            &old_face_indices.iter().map(|&i| i as u64).collect::<Vec<_>>(),
            new_face_idx as u64,
            rcad_kernel::persistent_naming::NamePropagationPolicy::Preserve,
        )
    }

    /// Get the number of recorded snapshots.
    pub fn snapshot_count(&self) -> usize {
        self.snapshots.len()
    }

    /// Check if undo is available.
    pub fn can_undo(&self) -> bool {
        self.current_snapshot > 0
    }

    /// Check if redo is available.
    pub fn can_redo(&self) -> bool {
        self.current_snapshot < self.snapshots.len().saturating_sub(1)
    }

    /// Get the action name at the current snapshot.
    pub fn current_action(&self) -> Option<&str> {
        self.snapshots.get(self.current_snapshot).map(|s| s.action.as_str())
    }

    /// Export naming events as a serializable history.
    pub fn export_naming_history(&self) -> rcad_kernel::persistent_naming::NamingHistory {
        self.naming_engine.export_naming_history()
    }
}

/// Convert a TopoNode to a unique entity ID.
fn node_to_entity_id(node: TopoNode) -> u64 {
    // Encode kind and index into a single u64.
    // Kind uses high 8 bits, index uses low 56 bits.
    let kind_bits = match node.kind {
        NodeKind::Solid => 0u64,
        NodeKind::Shell => 1u64,
        NodeKind::Face => 2u64,
        NodeKind::Wire => 3u64,
        NodeKind::Edge => 4u64,
        NodeKind::Vertex => 5u64,
    };
    (kind_bits << 56) | (node.index as u64)
}

/// Convert an entity ID back to a TopoNode (requires known kind).
fn entity_id_to_node(entity_id: u64, kind: NodeKind) -> TopoNode {
    let index = (entity_id & 0x00FFFFFFFFFFFFFF) as usize;
    TopoNode { kind, index }
}

/// NamedGraph: A TopoGraph with integrated naming history.
///
/// This provides a convenient wrapper for applications that need
/// automatic naming tracking during graph mutations.
#[derive(Debug, Clone)]
pub struct NamedGraph {
    graph: TopoGraph,
    history: BRepGraphHistory,
}

impl NamedGraph {
    /// Create a new named graph from a BRep.
    pub fn from_brep(brep: &BRep) -> Self {
        let graph = TopoGraph::from_brep(brep);
        let mut history = BRepGraphHistory::new();

        // Assign persistent IDs to all nodes.
        for node in &graph.nodes {
            history.assign_persistent_id(*node);
        }

        Self { graph, history }
    }

    /// Get the underlying graph.
    pub fn graph(&self) -> &TopoGraph {
        &self.graph
    }

    /// Get mutable access to the graph.
    pub fn graph_mut(&mut self) -> &mut TopoGraph {
        &mut self.graph
    }

    /// Get the history.
    pub fn history(&self) -> &BRepGraphHistory {
        &self.history
    }

    /// Get mutable access to the history.
    pub fn history_mut(&mut self) -> &mut BRepGraphHistory {
        &mut self.history
    }

    /// Apply a mutation with automatic naming tracking.
    pub fn mutate_tracked<F>(
        &mut self,
        action: &str,
        operation_type: OperationType,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph, &mut BRepGraphHistory),
    {
        self.history.begin_operation(operation_type, Some(action.to_string()));

        let before_nodes = self.graph.nodes.clone();
        mutator(&mut self.graph, &mut self.history);

        let issues = self.graph.validate();
        if issues.is_empty() {
            let after_nodes = self.graph.nodes.clone();
            let entity_map: Vec<(TopoNode, Option<TopoNode>)> = before_nodes
                .iter()
                .filter_map(|old| {
                    // Find if this node still exists.
                    let still_exists = after_nodes.contains(old);
                    if still_exists {
                        Some((*old, Some(*old)))
                    } else {
                        Some((*old, None))
                    }
                })
                .collect();

            self.history.record_mutation(&self.graph, action, &entity_map);

            self.history.finalize_operation(OperationStats {
                entity_count_before: before_nodes.len(),
                entity_count_after: after_nodes.len(),
                names_preserved: entity_map.iter().filter(|(_, new)| new.is_some()).count(),
                names_lost: entity_map.iter().filter(|(_, new)| new.is_none()).count(),
                names_generated: 0,
                conflicts_resolved: 0,
            });

            Ok(())
        } else {
            self.history.cancel_operation();
            Err(issues)
        }
    }

    /// Get the persistent ID for a node.
    pub fn get_persistent_id(&self, node: TopoNode) -> Option<PersistentId> {
        self.history.resolve_persistent(node)
    }

    /// Generate a stability report.
    pub fn stability_report(&self) -> CrossOperationStabilityReport {
        self.history.stability_report()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cross-Operation Naming Stability Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Metrics measuring naming stability for a single operation.
#[derive(Debug, Clone)]
pub struct OperationStabilityMetrics {
    /// Operation identifier.
    pub operation_id: rcad_kernel::persistent_naming::OperationId,
    /// Type of operation performed.
    pub operation_type: rcad_kernel::persistent_naming::OperationType,
    /// Optional label for the operation.
    pub label: Option<String>,
    /// Number of entities that retained their names through this operation.
    pub names_retained: usize,
    /// Number of entities that lost their names during this operation.
    pub names_lost: usize,
    /// Number of new names generated during this operation.
    pub names_generated: usize,
    /// Number of naming conflicts that occurred during this operation.
    pub conflicts: usize,
    /// Stability score for this specific operation (0.0 - 1.0).
    pub stability_score: f64,
    /// Cumulative stability score up to and including this operation.
    pub cumulative_stability: f64,
}

impl Default for OperationStabilityMetrics {
    fn default() -> Self {
        Self {
            operation_id: rcad_kernel::persistent_naming::OperationId::NULL,
            operation_type: rcad_kernel::persistent_naming::OperationType::Generic,
            label: None,
            names_retained: 0,
            names_lost: 0,
            names_generated: 0,
            conflicts: 0,
            stability_score: 1.0,
            cumulative_stability: 1.0,
        }
    }
}

/// Information about a broken naming chain.
#[derive(Debug, Clone)]
pub struct BrokenChainInfo {
    /// The persistent ID whose chain was broken.
    pub persistent_id: PersistentId,
    /// The operation where the break occurred.
    pub broken_at_operation: rcad_kernel::persistent_naming::OperationId,
    /// Entity ID that lost the name.
    pub entity_id: u64,
    /// Entity type if known.
    pub entity_type: Option<rcad_kernel::persistent_naming::EntityType>,
    /// Number of operations the chain survived before breaking.
    pub survived_operations: usize,
    /// Description of how the chain broke.
    pub break_reason: String,
}

/// A naming conflict detected across operations.
#[derive(Debug, Clone)]
pub struct NamingConflict {
    /// The persistent ID involved in the conflict.
    pub persistent_id: PersistentId,
    /// Operations where the conflict manifested.
    pub involved_operations: Vec<rcad_kernel::persistent_naming::OperationId>,
    /// Entity IDs that were in conflict.
    pub conflicting_entities: Vec<u64>,
    /// Type of conflict.
    pub conflict_type: ConflictType,
    /// Severity of the conflict.
    pub severity: rcad_kernel::persistent_naming::IssueSeverity,
    /// Whether the conflict was automatically resolved.
    pub auto_resolved: bool,
}

/// Types of naming conflicts that can occur.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConflictType {
    /// Same persistent ID assigned to multiple entities.
    DuplicateAssignment,
    /// Entity references a deleted persistent ID.
    ReferenceToDeleted,
    /// Genealogy chain is incomplete or broken.
    BrokenGenealogy,
    /// Unexpected name change during propagation.
    UnexpectedNameChange,
    /// Merge operation lost entity tracking.
    MergeTrackingLoss,
}

/// A recommendation for improving naming stability.
#[derive(Debug, Clone)]
pub struct StabilityRecommendation {
    /// Priority of this recommendation (higher = more important).
    pub priority: u32,
    /// Category of the recommendation.
    pub category: RecommendationCategory,
    /// Human-readable description of the recommendation.
    pub description: String,
    /// Operations this recommendation applies to (empty = all).
    pub affected_operations: Vec<rcad_kernel::persistent_naming::OperationId>,
    /// Estimated impact on stability score if implemented.
    pub estimated_impact: f64,
    /// Code or configuration suggestion.
    pub suggestion: Option<String>,
}

/// Categories of stability recommendations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecommendationCategory {
    /// Adjust naming rule selection.
    NamingRule,
    /// Improve name propagation policy.
    PropagationPolicy,
    /// Fix specific operation handling.
    OperationHandling,
    /// Improve entity tracking.
    EntityTracking,
    /// Address conflict resolution.
    ConflictResolution,
    /// General architecture improvement.
    Architecture,
}

/// Comprehensive analysis of naming stability across multiple operations.
#[derive(Debug, Clone)]
pub struct CrossOperationNamingAnalysis {
    /// Operations analyzed.
    pub operations: Vec<OperationRecord>,
    /// Entities tracked through all operations.
    pub entity_genealogy: HashMap<PersistentId, EntityGenealogy>,
    /// Stability metrics per operation.
    pub per_operation_stability: Vec<OperationStabilityMetrics>,
    /// Overall stability score (0.0 - 1.0).
    pub overall_stability: f64,
    /// Entities with broken naming chains.
    pub broken_chains: Vec<BrokenChainInfo>,
    /// Trend direction: positive = improving, negative = degrading.
    pub stability_trend: f64,
    /// Number of entities tracked at each operation boundary.
    pub entity_counts: Vec<usize>,
}

impl CrossOperationNamingAnalysis {
    /// Returns true if overall stability is excellent (> 95%).
    pub fn is_excellent(&self) -> bool {
        self.overall_stability >= 0.95
    }

    /// Returns true if overall stability is good (> 90%).
    pub fn is_good(&self) -> bool {
        self.overall_stability >= 0.90
    }

    /// Returns true if there are significant stability issues.
    pub fn has_issues(&self) -> bool {
        self.overall_stability < 0.90 || !self.broken_chains.is_empty()
    }

    /// Get the most problematic operation (lowest stability score).
    pub fn most_problematic_operation(&self) -> Option<&OperationStabilityMetrics> {
        self.per_operation_stability
            .iter()
            .min_by(|a, b| a.stability_score.partial_cmp(&b.stability_score).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Get operations sorted by stability score (ascending).
    pub fn operations_by_stability(&self) -> Vec<&OperationStabilityMetrics> {
        let mut ops: Vec<_> = self.per_operation_stability.iter().collect();
        ops.sort_by(|a, b| a.stability_score.partial_cmp(&b.stability_score).unwrap_or(std::cmp::Ordering::Equal));
        ops
    }

    /// Calculate the average stability score across all operations.
    pub fn average_operation_stability(&self) -> f64 {
        if self.per_operation_stability.is_empty() {
            return 1.0;
        }
        let sum: f64 = self.per_operation_stability.iter().map(|m| m.stability_score).sum();
        sum / self.per_operation_stability.len() as f64
    }

    /// Generate a human-readable summary.
    pub fn summary(&self) -> String {
        format!(
            "Cross-Operation Naming Analysis:\n\
             - Operations: {}\n\
             - Entities Tracked: {}\n\
             - Overall Stability: {:.1}%\n\
             - Stability Trend: {}\n\
             - Broken Chains: {}\n\
             - Avg Operation Stability: {:.1}%",
            self.operations.len(),
            self.entity_genealogy.len(),
            self.overall_stability * 100.0,
            if self.stability_trend > 0.0 { "Improving" } else if self.stability_trend < 0.0 { "Degrading" } else { "Stable" },
            self.broken_chains.len(),
            self.average_operation_stability() * 100.0
        )
    }
}

/// Track naming through a sequence of operations.
///
/// This function analyzes a series of BRepGraphHistory snapshots to determine
/// how naming stability evolves across multiple operations.
pub fn analyze_naming_sequence(
    history: &[BRepGraphHistory],
    initial_entities: &[TopoNode],
) -> CrossOperationNamingAnalysis {
    use rcad_kernel::persistent_naming::NamingEvent;

    if history.is_empty() {
        return CrossOperationNamingAnalysis {
            operations: Vec::new(),
            entity_genealogy: HashMap::new(),
            per_operation_stability: Vec::new(),
            overall_stability: 1.0,
            broken_chains: Vec::new(),
            stability_trend: 0.0,
            entity_counts: vec![initial_entities.len()],
        };
    }

    let mut operations: Vec<OperationRecord> = Vec::new();
    let mut entity_genealogy: HashMap<PersistentId, EntityGenealogy> = HashMap::new();
    let mut per_operation_stability: Vec<OperationStabilityMetrics> = Vec::new();
    let mut broken_chains: Vec<BrokenChainInfo> = Vec::new();
    let mut entity_counts: Vec<usize> = vec![initial_entities.len()];

    // Track entity persistence across operations.
    let mut entity_to_pid: HashMap<u64, PersistentId> = HashMap::new();
    let mut pid_to_entity: HashMap<PersistentId, u64> = HashMap::new();

    // Assign initial persistent IDs.
    for node in initial_entities {
        let entity_id = node_to_entity_id(*node);
        let pid = PersistentId(entity_id); // Use entity_id as basis for PID
        entity_to_pid.insert(entity_id, pid);
        pid_to_entity.insert(pid, entity_id);

        entity_genealogy.insert(pid, EntityGenealogy {
            persistent_id: pid,
            created_in_operation: rcad_kernel::persistent_naming::OperationId::NULL,
            evolution: vec![(rcad_kernel::persistent_naming::OperationId::NULL, entity_id)],
            current_entity_id: Some(entity_id),
            is_deleted: false,
        });
    }

    let mut cumulative_stability = 1.0;

    for hist in history.iter() {
        let cross_op = hist.naming_engine().cross_operation_history();

        // Copy operation records.
        for op in &cross_op.operations {
            operations.push(op.clone());
        }

        // Calculate stability metrics for each operation in this history.
        for op in &cross_op.operations {
            let total_entities = op.stats.entity_count_before.max(1);
            let preserved = op.stats.names_preserved;
            let stability_score = preserved as f64 / total_entities as f64;

            cumulative_stability *= stability_score;

            per_operation_stability.push(OperationStabilityMetrics {
                operation_id: op.id,
                operation_type: op.operation_type,
                label: op.label.clone(),
                names_retained: op.stats.names_preserved,
                names_lost: op.stats.names_lost,
                names_generated: op.stats.names_generated,
                conflicts: op.stats.conflicts_resolved,
                stability_score,
                cumulative_stability,
            });

            // Track entity counts.
            entity_counts.push(op.stats.entity_count_after);

            // Detect broken chains.
            for event in &op.naming_events {
                if let NamingEvent::Lost { entity_id, persistent_id } = event {
                    // Find how many operations this entity survived.
                    let survived_operations = entity_genealogy
                        .get(persistent_id)
                        .map(|g| g.evolution.len())
                        .unwrap_or(0);

                    broken_chains.push(BrokenChainInfo {
                        persistent_id: *persistent_id,
                        broken_at_operation: op.id,
                        entity_id: *entity_id,
                        entity_type: infer_entity_type_from_id(*entity_id),
                        survived_operations,
                        break_reason: "Entity removed without successor".to_string(),
                    });
                }
            }
        }

        // Update genealogy from this history.
        for (pid, genealogy) in cross_op.genealogy.iter() {
            entity_genealogy.insert(*pid, genealogy.clone());
        }
    }

    // Calculate overall stability.
    let total_preserved: usize = per_operation_stability.iter().map(|m| m.names_retained).sum();
    let total_lost: usize = per_operation_stability.iter().map(|m| m.names_lost).sum();
    let total = total_preserved + total_lost;
    let overall_stability = if total > 0 {
        total_preserved as f64 / total as f64
    } else {
        1.0
    };

    // Calculate stability trend.
    let stability_trend = calculate_stability_trend(&per_operation_stability);

    CrossOperationNamingAnalysis {
        operations,
        entity_genealogy,
        per_operation_stability,
        overall_stability,
        broken_chains,
        stability_trend,
        entity_counts,
    }
}

/// Detect naming conflicts across operations.
///
/// This function analyzes cross-operation data to identify conflicts that
/// may not be visible in single-operation analysis.
pub fn detect_cross_operation_conflicts(
    analysis: &CrossOperationNamingAnalysis,
) -> Vec<NamingConflict> {
    let mut conflicts: Vec<NamingConflict> = Vec::new();
    let mut pid_to_entities: HashMap<PersistentId, Vec<(u64, rcad_kernel::persistent_naming::OperationId)>> = HashMap::new();

    // Build a map of persistent ID to all entities that have held it.
    for (pid, genealogy) in &analysis.entity_genealogy {
        let mut entities = Vec::new();
        for (op_id, entity_id) in &genealogy.evolution {
            entities.push((*entity_id, *op_id));
        }
        pid_to_entities.insert(*pid, entities);
    }

    // Detect duplicate assignments (same PID to different entities at same time).
    for (pid, genealogy) in &analysis.entity_genealogy {
        if genealogy.is_deleted {
            continue;
        }

        // Check for entities that reference this PID after it was marked deleted.
        if let Some(current_entity) = genealogy.current_entity_id {
            // Verify the current entity still has this PID.
            for (other_pid, other_genealogy) in &analysis.entity_genealogy {
                if other_pid != pid && other_genealogy.current_entity_id == Some(current_entity) {
                    // Same entity has multiple PIDs.
                    conflicts.push(NamingConflict {
                        persistent_id: *pid,
                        involved_operations: vec![genealogy.created_in_operation],
                        conflicting_entities: vec![current_entity],
                        conflict_type: ConflictType::DuplicateAssignment,
                        severity: rcad_kernel::persistent_naming::IssueSeverity::Severe,
                        auto_resolved: false,
                    });
                }
            }
        }
    }

    // Detect broken genealogies.
    for (pid, genealogy) in &analysis.entity_genealogy {
        if genealogy.evolution.is_empty() && !genealogy.is_deleted {
            conflicts.push(NamingConflict {
                persistent_id: *pid,
                involved_operations: vec![genealogy.created_in_operation],
                conflicting_entities: vec![],
                conflict_type: ConflictType::BrokenGenealogy,
                severity: rcad_kernel::persistent_naming::IssueSeverity::Moderate,
                auto_resolved: false,
            });
        }
    }

    // Detect reference to deleted entities.
    for chain in &analysis.broken_chains {
        // Check if any genealogy still references this broken chain.
        if let Some(genealogy) = analysis.entity_genealogy.get(&chain.persistent_id)
            && !genealogy.is_deleted && genealogy.evolution.len() > chain.survived_operations {
                conflicts.push(NamingConflict {
                    persistent_id: chain.persistent_id,
                    involved_operations: vec![chain.broken_at_operation],
                    conflicting_entities: vec![chain.entity_id],
                    conflict_type: ConflictType::ReferenceToDeleted,
                    severity: rcad_kernel::persistent_naming::IssueSeverity::Critical,
                    auto_resolved: false,
                });
            }
    }

    conflicts
}

/// Generate recommendations for improving naming stability.
///
/// Based on the analysis, this function produces actionable recommendations
/// for improving the naming stability of future operations.
pub fn generate_stability_recommendations(
    analysis: &CrossOperationNamingAnalysis,
) -> Vec<StabilityRecommendation> {
    let mut recommendations: Vec<StabilityRecommendation> = Vec::new();

    // Check overall stability.
    if analysis.overall_stability < 0.5 {
        recommendations.push(StabilityRecommendation {
            priority: 100,
            category: RecommendationCategory::Architecture,
            description: "Critical naming stability issues detected. Consider reviewing the entire naming strategy.".to_string(),
            affected_operations: vec![],
            estimated_impact: 0.3,
            suggestion: Some("Enable Hybrid naming rule and Preserve propagation policy".to_string()),
        });
    } else if analysis.overall_stability < 0.8 {
        recommendations.push(StabilityRecommendation {
            priority: 80,
            category: RecommendationCategory::NamingRule,
            description: "Moderate naming stability degradation. Review naming rule configuration.".to_string(),
            affected_operations: vec![],
            estimated_impact: 0.15,
            suggestion: Some("Consider switching to HistoryTracking naming rule for better traceability".to_string()),
        });
    }

    // Check for problematic operations.
    if let Some(problematic) = analysis.most_problematic_operation()
        && problematic.stability_score < 0.7 {
            recommendations.push(StabilityRecommendation {
                priority: 90,
                category: RecommendationCategory::OperationHandling,
                description: format!(
                    "Operation {:?} has low stability score ({:.1}%). Review operation-specific handling.",
                    problematic.operation_type,
                    problematic.stability_score * 100.0
                ),
                affected_operations: vec![problematic.operation_id],
                estimated_impact: 0.2,
                suggestion: Some("Ensure entity mapping is correctly tracked during this operation".to_string()),
            });
        }

    // Check for broken chains.
    if !analysis.broken_chains.is_empty() {
        let severe_breaks: Vec<_> = analysis.broken_chains.iter()
            .filter(|c| c.survived_operations > 5)
            .collect();

        if !severe_breaks.is_empty() {
            recommendations.push(StabilityRecommendation {
                priority: 85,
                category: RecommendationCategory::EntityTracking,
                description: format!(
                    "{} long-lived entities lost their naming chains. Improve entity tracking during mutations.",
                    severe_breaks.len()
                ),
                affected_operations: severe_breaks.iter().map(|c| c.broken_at_operation).collect(),
                estimated_impact: 0.25,
                suggestion: Some("Implement explicit entity mapping during split/merge operations".to_string()),
            });
        }
    }

    // Check stability trend.
    if analysis.stability_trend < -0.1 {
        recommendations.push(StabilityRecommendation {
            priority: 75,
            category: RecommendationCategory::PropagationPolicy,
            description: "Stability is degrading over time. Naming propagation may need adjustment.".to_string(),
            affected_operations: vec![],
            estimated_impact: 0.1,
            suggestion: Some("Review NamePropagationPolicy settings for recent operations".to_string()),
        });
    }

    // Check for conflicts.
    let conflicts = detect_cross_operation_conflicts(analysis);
    let critical_conflicts: Vec<_> = conflicts.iter()
        .filter(|c| c.severity == rcad_kernel::persistent_naming::IssueSeverity::Critical)
        .collect();

    if !critical_conflicts.is_empty() {
        recommendations.push(StabilityRecommendation {
            priority: 95,
            category: RecommendationCategory::ConflictResolution,
            description: format!(
                "{} critical naming conflicts detected. Immediate resolution required.",
                critical_conflicts.len()
            ),
            affected_operations: critical_conflicts.iter().flat_map(|c| c.involved_operations.clone()).collect(),
            estimated_impact: 0.3,
            suggestion: Some("Manually resolve conflicts or reset naming for affected entities".to_string()),
        });
    }

    // Check per-operation patterns.
    let low_stability_ops: Vec<_> = analysis.per_operation_stability.iter()
        .filter(|m| m.stability_score < 0.8)
        .collect();

    if low_stability_ops.len() > analysis.per_operation_stability.len() / 2 {
        recommendations.push(StabilityRecommendation {
            priority: 70,
            category: RecommendationCategory::Architecture,
            description: "Multiple operations have low stability. Consider system-wide naming improvements.".to_string(),
            affected_operations: low_stability_ops.iter().map(|m| m.operation_id).collect(),
            estimated_impact: 0.2,
            suggestion: Some("Enable comprehensive entity tracking and increase name propagation fidelity".to_string()),
        });
    }

    // Sort by priority (highest first).
    recommendations.sort_by(|a, b| b.priority.cmp(&a.priority));
    recommendations
}

/// Calculate the stability trend from per-operation metrics.
///
/// Returns a value between -1.0 (strongly degrading) and 1.0 (strongly improving).
fn calculate_stability_trend(metrics: &[OperationStabilityMetrics]) -> f64 {
    if metrics.len() < 2 {
        return 0.0;
    }

    // Simple linear regression on stability scores.
    let n = metrics.len() as f64;
    let sum_x: f64 = (0..metrics.len()).map(|i| i as f64).sum();
    let sum_y: f64 = metrics.iter().map(|m| m.stability_score).sum();
    let sum_xy: f64 = metrics.iter().enumerate()
        .map(|(i, m)| i as f64 * m.stability_score)
        .sum();
    let sum_x2: f64 = (0..metrics.len()).map(|i| (i * i) as f64).sum();

    let denominator = n * sum_x2 - sum_x * sum_x;
    if denominator.abs() < f64::EPSILON {
        return 0.0;
    }

    let slope = (n * sum_xy - sum_x * sum_y) / denominator;

    // Normalize slope to [-1, 1] range.
    // A slope of 0.1 per operation is considered a strong trend.
    slope.clamp(-1.0, 1.0)
}

/// Infer entity type from an encoded entity ID.
fn infer_entity_type_from_id(entity_id: u64) -> Option<rcad_kernel::persistent_naming::EntityType> {
    // High 8 bits encode the kind.
    let kind_bits = entity_id >> 56;
    match kind_bits {
        0 | 1 => Some(rcad_kernel::persistent_naming::EntityType::Solid),
        2 => Some(rcad_kernel::persistent_naming::EntityType::Face),
        4 => Some(rcad_kernel::persistent_naming::EntityType::Edge),
        5 => Some(rcad_kernel::persistent_naming::EntityType::Vertex),
        _ => None,
    }
}

// Re-export types needed for analysis.
pub use rcad_kernel::persistent_naming::{
    OperationRecord, EntityGenealogy,
};

impl BRepGraphHistory {
    /// Cancel the current operation.
    fn cancel_operation(&mut self) {
        self.naming_engine.cancel_operation();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced Persistent Naming Semantics
// ─────────────────────────────────────────────────────────────────────────────

/// A scoped identifier for a topological entity within a naming context.
///
/// `ScopedId` combines a persistent ID with a naming scope (part, assembly, operation)
/// to provide fully-qualified identifiers that are unique across an entire model hierarchy.
///
/// # Example
/// ```
/// use rcad_algorithms::brep_graph::{ScopedId, NamingScope};
/// use rcad_kernel::PersistentId;
///
/// // A face with ID 42 in part "housing", assembly "machine", operation "fillet"
/// let scoped = ScopedId {
///     persistent_id: PersistentId(42),
///     scope: NamingScope {
///         part: Some("housing".to_string()),
///         assembly: Some("machine".to_string()),
///         operation: Some("fillet".to_string()),
///     },
/// };
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ScopedId {
    /// The stable persistent ID for this entity.
    pub persistent_id: PersistentId,
    /// The naming scope in which this ID is defined.
    pub scope: NamingScope,
}

impl ScopedId {
    /// Create a new scoped ID with the given persistent ID and scope.
    pub fn new(persistent_id: PersistentId, scope: NamingScope) -> Self {
        Self { persistent_id, scope }
    }

    /// Create a scoped ID with a null persistent ID and empty scope.
    pub fn null() -> Self {
        Self {
            persistent_id: PersistentId::NULL,
            scope: NamingScope::default(),
        }
    }

    /// Returns true if this is a null/invalid scoped ID.
    pub fn is_null(&self) -> bool {
        self.persistent_id.is_null()
    }

    /// Generate a fully-qualified name string for this scoped ID.
    pub fn qualified_name(&self) -> String {
        let mut parts = Vec::new();
        if let Some(ref assembly) = self.scope.assembly {
            parts.push(assembly.clone());
        }
        if let Some(ref part) = self.scope.part {
            parts.push(part.clone());
        }
        if let Some(ref op) = self.scope.operation {
            parts.push(op.clone());
        }
        parts.push(format!("e{}", self.persistent_id.raw()));
        parts.join("::")
    }
}

/// The naming scope defines the context in which persistent IDs are meaningful.
///
/// Scopes form a hierarchy: assembly > part > operation. An entity's full identity
/// is determined by its persistent ID within the current scope.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NamingScope {
    /// The part name (e.g., "housing", "cover").
    pub part: Option<String>,
    /// The assembly name (e.g., "machine", "device").
    pub assembly: Option<String>,
    /// The operation that created or last modified this entity.
    pub operation: Option<String>,
}

impl NamingScope {
    /// Create a new empty scope.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a scope for a specific part.
    pub fn for_part(part: impl Into<String>) -> Self {
        Self {
            part: Some(part.into()),
            assembly: None,
            operation: None,
        }
    }

    /// Create a scope for a specific assembly.
    pub fn for_assembly(assembly: impl Into<String>) -> Self {
        Self {
            part: None,
            assembly: Some(assembly.into()),
            operation: None,
        }
    }

    /// Create a scope for a specific operation within a part.
    pub fn for_operation(part: impl Into<String>, operation: impl Into<String>) -> Self {
        Self {
            part: Some(part.into()),
            assembly: None,
            operation: Some(operation.into()),
        }
    }

    /// Set the part name.
    pub fn with_part(mut self, part: impl Into<String>) -> Self {
        self.part = Some(part.into());
        self
    }

    /// Set the assembly name.
    pub fn with_assembly(mut self, assembly: impl Into<String>) -> Self {
        self.assembly = Some(assembly.into());
        self
    }

    /// Set the operation name.
    pub fn with_operation(mut self, operation: impl Into<String>) -> Self {
        self.operation = Some(operation.into());
        self
    }

    /// Create a child scope for a sub-operation.
    pub fn child_scope(&self, operation: impl Into<String>) -> Self {
        Self {
            part: self.part.clone(),
            assembly: self.assembly.clone(),
            operation: Some(operation.into()),
        }
    }

    /// Check if this scope is a parent of (or equal to) another scope.
    pub fn contains(&self, other: &NamingScope) -> bool {
        match (&self.assembly, &other.assembly) {
            (Some(a1), Some(a2)) if a1 != a2 => return false,
            (Some(_), None) => return false,
            _ => {}
        }
        match (&self.part, &other.part) {
            (Some(p1), Some(p2)) if p1 != p2 => return false,
            (Some(_), None) => return false,
            _ => {}
        }
        true
    }
}

/// Enhanced naming context that tracks scopes and entity relationships.
///
/// `EnhancedNamingContext` extends the basic `NamingContext` with:
/// - Scope-aware ID assignment
/// - Detailed genealogy tracking
/// - Conflict detection and resolution
/// - Serialization support
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EnhancedNamingContext {
    /// The current naming scope.
    pub current_scope: NamingScope,
    /// Mapping from entity IDs to scoped IDs.
    entity_to_scoped: HashMap<u64, ScopedId>,
    /// Reverse mapping from scoped IDs to entity IDs.
    scoped_to_entity: HashMap<ScopedId, u64>,
    /// Genealogy records indexed by persistent ID.
    genealogy: HashMap<PersistentId, EntityGenealogyRecord>,
    /// Pending name assignments waiting for scope resolution.
    pending_assignments: Vec<PendingNameAssignment>,
    /// Conflict resolution history.
    conflict_history: Vec<NameConflictRecord>,
    /// Next persistent ID to allocate.
    next_persistent_id: u64,
}

/// Record of an entity's genealogy through operations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntityGenealogyRecord {
    /// The persistent ID being tracked.
    pub persistent_id: PersistentId,
    /// The scope in which this entity was created.
    pub creation_scope: NamingScope,
    /// The operation that created this entity.
    pub creation_operation: Option<String>,
    /// Chain of transformations: (operation, entity_id_before, entity_id_after).
    pub transformation_chain: Vec<GenealogyStep>,
    /// Parent entity IDs (for merged entities, this has multiple entries).
    pub parent_ids: Vec<PersistentId>,
    /// Child entity IDs (for split entities, this has multiple entries).
    pub child_ids: Vec<PersistentId>,
    /// Current status of this entity.
    pub status: EntityStatus,
}

/// A single step in an entity's genealogy transformation chain.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenealogyStep {
    /// The operation that caused this transformation.
    pub operation: String,
    /// The entity ID before this operation (None if generated).
    pub entity_id_before: Option<u64>,
    /// The entity ID after this operation.
    pub entity_id_after: u64,
    /// The scope at the time of this operation.
    pub scope: NamingScope,
}

/// Status of an entity in the genealogy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EntityStatus {
    /// Entity is active and present in the model.
    Active,
    /// Entity was deleted or consumed by an operation.
    Deleted,
    /// Entity was merged into another entity.
    Merged,
    /// Entity was split into multiple entities.
    Split,
    /// Entity is pending resolution.
    Pending,
}

/// A pending name assignment waiting for scope resolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PendingNameAssignment {
    /// The entity ID to be assigned.
    pub entity_id: u64,
    /// The proposed scope for this assignment.
    pub proposed_scope: NamingScope,
    /// Source entity IDs this entity was derived from.
    pub source_entities: Vec<u64>,
    /// The propagation policy to use.
    pub propagation_policy: NamePropagationPolicy,
}

/// Record of a naming conflict and its resolution.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NameConflictRecord {
    /// The conflicting persistent ID.
    pub persistent_id: PersistentId,
    /// The entity IDs involved in the conflict.
    pub conflicting_entities: Vec<u64>,
    /// The operation where the conflict occurred.
    pub operation: String,
    /// The scope where the conflict occurred.
    pub scope: NamingScope,
    /// How the conflict was resolved.
    pub resolution: NameConflictResolution,
    /// Timestamp of the conflict (sequence number).
    pub sequence: u64,
}

/// Strategies for resolving naming conflicts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NameConflictResolution {
    /// Kept the existing binding, rejected the new.
    KeepExisting,
    /// Replaced the existing binding with the new.
    ReplaceWithNew,
    /// Generated a new persistent ID for the new entity.
    GenerateNewId,
    /// Merged both entities under a shared context.
    MergeEntities,
    /// Created an alias mapping.
    CreateAlias,
    /// Could not resolve automatically - requires manual intervention.
    Unresolved,
}

/// Operation-specific name propagation rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NamePropagationRule {
    /// The operation type this rule applies to.
    pub operation_type: OperationType,
    /// Policy for face entities.
    pub face_policy: NamePropagationPolicy,
    /// Policy for edge entities.
    pub edge_policy: NamePropagationPolicy,
    /// Policy for vertex entities.
    pub vertex_policy: NamePropagationPolicy,
    /// Whether to track genealogy for this operation.
    pub track_genealogy: bool,
    /// Conflict resolution strategy for this operation.
    pub conflict_resolution: NameConflictResolution,
}

impl NamePropagationRule {
    /// Create a default propagation rule for an operation type.
    pub fn for_operation(operation_type: OperationType) -> Self {
        match operation_type {
            OperationType::BooleanUnion |
            OperationType::BooleanIntersection |
            OperationType::BooleanDifference => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Preserve,
                edge_policy: NamePropagationPolicy::Preserve,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::GenerateNewId,
            },
            OperationType::Feature => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Inherit,
                edge_policy: NamePropagationPolicy::Inherit,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::KeepExisting,
            },
            OperationType::EdgeSplit |
            OperationType::FaceSplit => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Inherit,
                edge_policy: NamePropagationPolicy::Inherit,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::GenerateNewId,
            },
            OperationType::Merge => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Combine,
                edge_policy: NamePropagationPolicy::Combine,
                vertex_policy: NamePropagationPolicy::Combine,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::MergeEntities,
            },
            OperationType::Delete => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Generate,
                edge_policy: NamePropagationPolicy::Generate,
                vertex_policy: NamePropagationPolicy::Generate,
                track_genealogy: false,
                conflict_resolution: NameConflictResolution::KeepExisting,
            },
            OperationType::Transform => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Preserve,
                edge_policy: NamePropagationPolicy::Preserve,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::KeepExisting,
            },
            OperationType::Generic |
            OperationType::Import => Self {
                operation_type,
                face_policy: NamePropagationPolicy::Preserve,
                edge_policy: NamePropagationPolicy::Preserve,
                vertex_policy: NamePropagationPolicy::Preserve,
                track_genealogy: true,
                conflict_resolution: NameConflictResolution::GenerateNewId,
            },
        }
    }

    /// Get the propagation policy for a specific entity kind.
    pub fn policy_for_kind(&self, kind: NodeKind) -> NamePropagationPolicy {
        match kind {
            NodeKind::Face | NodeKind::Shell | NodeKind::Solid => self.face_policy,
            NodeKind::Edge | NodeKind::Wire => self.edge_policy,
            NodeKind::Vertex => self.vertex_policy,
        }
    }
}

/// Propagation policies for name inheritance through operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[derive(Default)]
pub enum NamePropagationPolicy {
    /// Keep the original entity's name unchanged.
    #[default]
    Preserve,
    /// Inherit the parent entity's name with a disambiguating suffix.
    Inherit,
    /// Generate a completely new name.
    Generate,
    /// Combine names from multiple source entities (for merges).
    Combine,
    /// Create a derivative name based on geometric properties.
    GeometryBased,
    /// Create a derivative name based on topological relationships.
    TopologyBased,
}
include!("e1.rs");

