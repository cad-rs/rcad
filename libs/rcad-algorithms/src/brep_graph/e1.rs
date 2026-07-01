

impl EnhancedNamingContext {
    /// Create a new empty naming context.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a naming context with a specific scope.
    pub fn with_scope(scope: NamingScope) -> Self {
        Self {
            current_scope: scope,
            ..Default::default()
        }
    }

    /// Set the current naming scope.
    pub fn set_scope(&mut self, scope: NamingScope) {
        self.current_scope = scope;
    }

    /// Get the current naming scope.
    pub fn scope(&self) -> &NamingScope {
        &self.current_scope
    }

    /// Assign a new persistent ID to an entity.
    pub fn assign_id(&mut self, entity_id: u64) -> PersistentId {
        let pid = self.allocate_persistent_id();
        let scoped_id = ScopedId::new(pid, self.current_scope.clone());
        self.entity_to_scoped.insert(entity_id, scoped_id.clone());
        self.scoped_to_entity.insert(scoped_id, entity_id);

        // Create genealogy record.
        self.genealogy.insert(pid, EntityGenealogyRecord {
            persistent_id: pid,
            creation_scope: self.current_scope.clone(),
            creation_operation: self.current_scope.operation.clone(),
            transformation_chain: vec![],
            parent_ids: vec![],
            child_ids: vec![],
            status: EntityStatus::Active,
        });

        pid
    }

    /// Assign a persistent ID derived from source entities.
    pub fn assign_derived_id(
        &mut self,
        entity_id: u64,
        source_entities: &[u64],
        policy: NamePropagationPolicy,
    ) -> PersistentId {
        match policy {
            NamePropagationPolicy::Preserve | NamePropagationPolicy::Inherit => {
                // Inherit from the first source that has a persistent ID.
                if let Some(&source_id) = source_entities.first()
                    && let Some(scoped) = self.entity_to_scoped.get(&source_id) {
                        let pid = scoped.persistent_id;
                        let new_scoped = ScopedId::new(pid, self.current_scope.clone());
                        self.entity_to_scoped.insert(entity_id, new_scoped.clone());
                        self.scoped_to_entity.insert(new_scoped, entity_id);

                        // Update genealogy.
                        if let Some(record) = self.genealogy.get_mut(&pid) {
                            record.transformation_chain.push(GenealogyStep {
                                operation: self.current_scope.operation.clone().unwrap_or_default(),
                                entity_id_before: Some(source_id),
                                entity_id_after: entity_id,
                                scope: self.current_scope.clone(),
                            });
                        }

                        return pid;
                    }
                self.assign_id(entity_id)
            }
            NamePropagationPolicy::Combine => {
                // Combine all source persistent IDs into the genealogy.
                let pid = self.assign_id(entity_id);
                if let Some(record) = self.genealogy.get_mut(&pid) {
                    for &source_id in source_entities {
                        if let Some(scoped) = self.entity_to_scoped.get(&source_id) {
                            record.parent_ids.push(scoped.persistent_id);
                        }
                    }
                }
                pid
            }
            NamePropagationPolicy::Generate |
            NamePropagationPolicy::GeometryBased |
            NamePropagationPolicy::TopologyBased => {
                self.assign_id(entity_id)
            }
        }
    }

    /// Resolve a persistent ID to an entity ID within the current scope.
    pub fn resolve_entity(&self, pid: PersistentId) -> Option<u64> {
        let scoped = ScopedId::new(pid, self.current_scope.clone());
        self.scoped_to_entity.get(&scoped).copied()
    }

    /// Resolve an entity ID to a persistent ID.
    pub fn resolve_persistent(&self, entity_id: u64) -> Option<PersistentId> {
        self.entity_to_scoped.get(&entity_id).map(|s| s.persistent_id)
    }

    /// Record a split operation: one entity becomes multiple.
    pub fn record_split(
        &mut self,
        source_entity_id: u64,
        target_entity_ids: &[u64],
        operation: &str,
    ) -> Vec<PersistentId> {
        let source_pid = self.resolve_persistent(source_entity_id);
        let mut result_pids = Vec::with_capacity(target_entity_ids.len());

        for (i, &target_id) in target_entity_ids.iter().enumerate() {
            let pid = if i == 0 {
                // First target inherits the source's persistent ID.
                if let Some(pid) = source_pid {
                    let scoped = ScopedId::new(pid, self.current_scope.clone());
                    self.entity_to_scoped.insert(target_id, scoped.clone());
                    self.scoped_to_entity.insert(scoped, target_id);
                    pid
                } else {
                    self.assign_id(target_id)
                }
            } else {
                // Subsequent targets get new IDs.
                self.assign_id(target_id)
            };
            result_pids.push(pid);
        }

        // Update genealogy for the source entity.
        if let Some(pid) = source_pid
            && let Some(record) = self.genealogy.get_mut(&pid) {
                record.status = EntityStatus::Split;
                record.child_ids.extend_from_slice(&result_pids[1..]);
                record.transformation_chain.push(GenealogyStep {
                    operation: operation.to_string(),
                    entity_id_before: Some(source_entity_id),
                    entity_id_after: result_pids.first().map(|&p| {
                        self.scoped_to_entity.get(&ScopedId::new(p, self.current_scope.clone()))
                            .copied()
                            .unwrap_or(0)
                    }).unwrap_or(0),
                    scope: self.current_scope.clone(),
                });
            }

        result_pids
    }

    /// Record a merge operation: multiple entities become one.
    pub fn record_merge(
        &mut self,
        source_entity_ids: &[u64],
        target_entity_id: u64,
        operation: &str,
        resolution: NameConflictResolution,
    ) -> PersistentId {
        // Find the first source with a persistent ID.
        let primary_pid = source_entity_ids
            .iter()
            .find_map(|&id| self.resolve_persistent(id));

        let target_pid = match resolution {
            NameConflictResolution::KeepExisting => {
                if let Some(pid) = primary_pid {
                    let scoped = ScopedId::new(pid, self.current_scope.clone());
                    self.entity_to_scoped.insert(target_entity_id, scoped.clone());
                    self.scoped_to_entity.insert(scoped, target_entity_id);
                    pid
                } else {
                    self.assign_id(target_entity_id)
                }
            }
            NameConflictResolution::GenerateNewId => {
                self.assign_id(target_entity_id)
            }
            NameConflictResolution::MergeEntities => {
                let pid = self.assign_id(target_entity_id);
                // Collect parent IDs first to avoid borrow conflict
                let parent_pids: Vec<PersistentId> = source_entity_ids
                    .iter()
                    .filter_map(|&source_id| self.resolve_persistent(source_id))
                    .collect();
                if let Some(record) = self.genealogy.get_mut(&pid) {
                    for source_pid in parent_pids {
                        record.parent_ids.push(source_pid);
                    }
                }
                pid
            }
            _ => {
                if let Some(pid) = primary_pid {
                    pid
                } else {
                    self.assign_id(target_entity_id)
                }
            }
        };

        // Mark source entities as merged.
        for &source_id in source_entity_ids {
            if let Some(pid) = self.resolve_persistent(source_id)
                && let Some(record) = self.genealogy.get_mut(&pid) {
                    record.status = EntityStatus::Merged;
                    record.child_ids.push(target_pid);
                    record.transformation_chain.push(GenealogyStep {
                        operation: operation.to_string(),
                        entity_id_before: Some(source_id),
                        entity_id_after: target_entity_id,
                        scope: self.current_scope.clone(),
                    });
                }
        }

        target_pid
    }

    /// Detect naming conflicts in the current state.
    pub fn detect_conflicts(&self) -> Vec<NameConflictRecord> {
        let mut conflicts = Vec::new();
        let mut pid_to_entities: HashMap<PersistentId, Vec<u64>> = HashMap::new();

        // Build a map of persistent ID to all entities that have it.
        for (&entity_id, scoped) in &self.entity_to_scoped {
            pid_to_entities
                .entry(scoped.persistent_id)
                .or_default()
                .push(entity_id);
        }

        // Find persistent IDs assigned to multiple active entities.
        for (pid, entities) in pid_to_entities {
            if entities.len() > 1 {
                // Check if all entities are active.
                let active_count = entities.iter()
                    .filter(|&&_entity_id| {
                        self.genealogy.get(&pid)
                            .map(|r| r.status == EntityStatus::Active)
                            .unwrap_or(false)
                    })
                    .count();

                if active_count > 1 {
                    conflicts.push(NameConflictRecord {
                        persistent_id: pid,
                        conflicting_entities: entities,
                        operation: self.current_scope.operation.clone().unwrap_or_default(),
                        scope: self.current_scope.clone(),
                        resolution: NameConflictResolution::Unresolved,
                        sequence: self.conflict_history.len() as u64,
                    });
                }
            }
        }

        conflicts
    }

    /// Resolve a naming conflict.
    pub fn resolve_conflict(
        &mut self,
        conflict: &NameConflictRecord,
        resolution: NameConflictResolution,
    ) -> Result<(), String> {
        match resolution {
            NameConflictResolution::KeepExisting => {
                // No action needed - first entity keeps the ID.
            }
            NameConflictResolution::ReplaceWithNew => {
                // Replace all but the last entity.
                for &entity_id in &conflict.conflicting_entities[..conflict.conflicting_entities.len() - 1] {
                    self.assign_id(entity_id);
                }
            }
            NameConflictResolution::GenerateNewId => {
                // Generate new IDs for all conflicting entities.
                for &entity_id in &conflict.conflicting_entities {
                    self.assign_id(entity_id);
                }
            }
            NameConflictResolution::MergeEntities => {
                // This is handled at the operation level.
            }
            NameConflictResolution::CreateAlias => {
                // Create alias mappings (not fully implemented here).
            }
            NameConflictResolution::Unresolved => {
                return Err("Conflict could not be resolved automatically".to_string());
            }
        }

        // Record the conflict resolution.
        let mut resolved = conflict.clone();
        resolved.resolution = resolution;
        self.conflict_history.push(resolved);

        Ok(())
    }

    /// Get the genealogy record for a persistent ID.
    pub fn get_genealogy(&self, pid: PersistentId) -> Option<&EntityGenealogyRecord> {
        self.genealogy.get(&pid)
    }

    /// Trace the full ancestry of an entity.
    pub fn trace_ancestry(&self, pid: PersistentId) -> Vec<PersistentId> {
        let mut ancestors = Vec::new();
        let mut visited = std::collections::HashSet::new();
        self.collect_ancestry(pid, &mut ancestors, &mut visited);
        ancestors
    }

    fn collect_ancestry(
        &self,
        pid: PersistentId,
        ancestors: &mut Vec<PersistentId>,
        visited: &mut std::collections::HashSet<PersistentId>,
    ) {
        if !visited.insert(pid) {
            return;
        }
        if let Some(record) = self.genealogy.get(&pid) {
            for &parent_pid in &record.parent_ids {
                ancestors.push(parent_pid);
                self.collect_ancestry(parent_pid, ancestors, visited);
            }
        }
    }

    /// Trace the full descendants of an entity.
    pub fn trace_descendants(&self, pid: PersistentId) -> Vec<PersistentId> {
        let mut descendants = Vec::new();
        let mut visited = std::collections::HashSet::new();
        self.collect_descendants(pid, &mut descendants, &mut visited);
        descendants
    }

    fn collect_descendants(
        &self,
        pid: PersistentId,
        descendants: &mut Vec<PersistentId>,
        visited: &mut std::collections::HashSet<PersistentId>,
    ) {
        if !visited.insert(pid) {
            return;
        }
        if let Some(record) = self.genealogy.get(&pid) {
            for &child_pid in &record.child_ids {
                descendants.push(child_pid);
                self.collect_descendants(child_pid, descendants, visited);
            }
        }
    }

    /// Mark an entity as deleted.
    pub fn mark_deleted(&mut self, entity_id: u64) {
        if let Some(pid) = self.resolve_persistent(entity_id)
            && let Some(record) = self.genealogy.get_mut(&pid) {
                record.status = EntityStatus::Deleted;
            }
    }

    /// Get all entities with a specific status.
    pub fn entities_by_status(&self, status: EntityStatus) -> Vec<PersistentId> {
        self.genealogy
            .iter()
            .filter(|(_, r)| r.status == status)
            .map(|(&pid, _)| pid)
            .collect()
    }

    /// Get the number of active entities.
    pub fn active_entity_count(&self) -> usize {
        self.entities_by_status(EntityStatus::Active).len()
    }

    /// Allocate a new persistent ID.
    fn allocate_persistent_id(&mut self) -> PersistentId {
        self.next_persistent_id += 1;
        PersistentId(self.next_persistent_id)
    }

    /// Clear all bindings and reset the context.
    pub fn clear(&mut self) {
        self.entity_to_scoped.clear();
        self.scoped_to_entity.clear();
        self.genealogy.clear();
        self.pending_assignments.clear();
        self.conflict_history.clear();
        self.next_persistent_id = 0;
    }

    /// Export the context state for serialization.
    pub fn export_state(&self) -> EnhancedNamingContextState {
        EnhancedNamingContextState {
            current_scope: self.current_scope.clone(),
            entity_to_scoped: self.entity_to_scoped.clone(),
            genealogy: self.genealogy.clone(),
            conflict_history: self.conflict_history.clone(),
            next_persistent_id: self.next_persistent_id,
        }
    }

    /// Import context state from a serialized form.
    pub fn import_state(&mut self, state: EnhancedNamingContextState) {
        self.current_scope = state.current_scope;
        self.entity_to_scoped = state.entity_to_scoped;
        self.genealogy = state.genealogy;
        self.conflict_history = state.conflict_history;
        self.next_persistent_id = state.next_persistent_id;

        // Rebuild reverse mapping.
        self.scoped_to_entity.clear();
        for (&entity_id, scoped) in &self.entity_to_scoped {
            self.scoped_to_entity.insert(scoped.clone(), entity_id);
        }
    }
}

/// Serializable state of an EnhancedNamingContext.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnhancedNamingContextState {
    pub current_scope: NamingScope,
    pub entity_to_scoped: HashMap<u64, ScopedId>,
    pub genealogy: HashMap<PersistentId, EntityGenealogyRecord>,
    pub conflict_history: Vec<NameConflictRecord>,
    pub next_persistent_id: u64,
}

/// Manager for operation-specific name propagation rules.
#[derive(Debug, Clone)]
pub struct NamePropagationManager {
    /// Rules indexed by operation type.
    rules: HashMap<OperationType, NamePropagationRule>,
    /// Default rule for unknown operation types.
    default_rule: NamePropagationRule,
}

impl Default for NamePropagationManager {
    fn default() -> Self {
        Self::new()
    }
}

impl NamePropagationManager {
    /// Create a new propagation manager with default rules.
    pub fn new() -> Self {
        let mut rules = HashMap::new();

        // Create default rules for all operation types.
        for op_type in [
            OperationType::BooleanUnion,
            OperationType::BooleanIntersection,
            OperationType::BooleanDifference,
            OperationType::EdgeSplit,
            OperationType::FaceSplit,
            OperationType::Merge,
            OperationType::Delete,
            OperationType::Transform,
            OperationType::Feature,
            OperationType::Generic,
            OperationType::Import,
        ] {
            rules.insert(op_type, NamePropagationRule::for_operation(op_type));
        }

        Self {
            rules,
            default_rule: NamePropagationRule::for_operation(OperationType::Generic),
        }
    }

    /// Get the propagation rule for an operation type.
    pub fn get_rule(&self, operation_type: OperationType) -> &NamePropagationRule {
        self.rules.get(&operation_type).unwrap_or(&self.default_rule)
    }

    /// Set a custom propagation rule for an operation type.
    pub fn set_rule(&mut self, rule: NamePropagationRule) {
        self.rules.insert(rule.operation_type, rule);
    }

    /// Apply a propagation rule to an entity transformation.
    pub fn apply_propagation(
        &self,
        context: &mut EnhancedNamingContext,
        operation_type: OperationType,
        source_entities: &[u64],
        target_entities: &[u64],
        entity_kind: NodeKind,
        operation_name: &str,
    ) -> Vec<PersistentId> {
        let rule = self.get_rule(operation_type);
        let policy = rule.policy_for_kind(entity_kind);

        // Handle split (1 -> many).
        if source_entities.len() == 1 && target_entities.len() > 1 {
            return context.record_split(source_entities[0], target_entities, operation_name);
        }

        // Handle merge (many -> 1).
        if source_entities.len() > 1 && target_entities.len() == 1 {
            return vec![context.record_merge(
                source_entities,
                target_entities[0],
                operation_name,
                rule.conflict_resolution,
            )];
        }

        // Handle 1 -> 1 transformation.
        if source_entities.len() == 1 && target_entities.len() == 1 {
            let pid = context.assign_derived_id(target_entities[0], source_entities, policy);
            return vec![pid];
        }

        // Handle generation (0 -> many).
        if source_entities.is_empty() {
            return target_entities.iter()
                .map(|&entity_id| context.assign_id(entity_id))
                .collect();
        }

        // Default: assign new IDs.
        target_entities.iter()
            .map(|&entity_id| context.assign_id(entity_id))
            .collect()
    }
}

/// Extension trait for BRepGraphHistory to support enhanced naming.
pub trait BRepGraphHistoryExt {
    /// Get the enhanced naming context.
    fn enhanced_context(&self) -> &EnhancedNamingContext;

    /// Get mutable access to the enhanced naming context.
    fn enhanced_context_mut(&mut self) -> &mut EnhancedNamingContext;

    /// Begin an operation with enhanced naming support.
    fn begin_enhanced_operation(
        &mut self,
        operation_type: OperationType,
        label: Option<String>,
        scope: NamingScope,
    );

    /// Propagate names through a boolean operation.
    fn propagate_boolean_names(
        &mut self,
        source_a_entities: &[TopoNode],
        source_b_entities: &[TopoNode],
        result_entities: &[TopoNode],
        operation: BooleanOperationType,
    );

    /// Propagate names through a fillet operation.
    fn propagate_fillet_names(
        &mut self,
        source_edges: &[usize],
        affected_faces: &[usize],
        new_faces: &[usize],
    );

    /// Propagate names through a chamfer operation.
    fn propagate_chamfer_names(
        &mut self,
        source_edges: &[usize],
        affected_faces: &[usize],
        new_faces: &[usize],
    );
}

/// Types of boolean operations for naming propagation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOperationType {
    Union,
    Intersection,
    Difference,
}

// ─────────────────────────────────────────────────────────────────────────────
// Serialization Support
// ─────────────────────────────────────────────────────────────────────────────

/// Serializable snapshot of a naming context for undo/redo support.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamingContextSnapshot {
    /// Unique identifier for this snapshot.
    pub id: u64,
    /// Timestamp when the snapshot was created.
    pub timestamp: u64,
    /// The scope at the time of the snapshot.
    pub scope: NamingScope,
    /// All entity-to-persistent-ID mappings.
    pub mappings: Vec<(u64, PersistentId, NamingScope)>,
    /// Genealogy records.
    pub genealogy: Vec<EntityGenealogyRecord>,
    /// Conflict records.
    pub conflicts: Vec<NameConflictRecord>,
    /// Operation that created this snapshot.
    pub operation: Option<String>,
}

impl NamingContextSnapshot {
    /// Create a snapshot from an enhanced naming context.
    pub fn from_context(context: &EnhancedNamingContext, id: u64, operation: Option<String>) -> Self {
        let mappings = context.entity_to_scoped.iter()
            .map(|(&entity_id, scoped)| {
                (entity_id, scoped.persistent_id, scoped.scope.clone())
            })
            .collect();

        Self {
            id,
            timestamp: std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0),
            scope: context.current_scope.clone(),
            mappings,
            genealogy: context.genealogy.values().cloned().collect(),
            conflicts: context.conflict_history.clone(),
            operation,
        }
    }

    /// Restore a naming context from this snapshot.
    pub fn restore_to(&self, context: &mut EnhancedNamingContext) {
        context.clear();
        context.current_scope = self.scope.clone();

        for (entity_id, pid, scope) in &self.mappings {
            let scoped = ScopedId::new(*pid, scope.clone());
            context.entity_to_scoped.insert(*entity_id, scoped.clone());
            context.scoped_to_entity.insert(scoped, *entity_id);
        }

        for record in &self.genealogy {
            context.genealogy.insert(record.persistent_id, record.clone());
        }

        context.conflict_history = self.conflicts.clone();
        context.next_persistent_id = self.genealogy.iter()
            .map(|r| r.persistent_id.raw())
            .max()
            .unwrap_or(0);
    }
}

/// Manager for naming context snapshots supporting undo/redo.
#[derive(Debug, Clone, Default)]
pub struct NamingSnapshotManager {
    /// All snapshots in chronological order.
    snapshots: Vec<NamingContextSnapshot>,
    /// Current position in the snapshot history.
    current_index: usize,
    /// Next snapshot ID.
    next_id: u64,
}

impl NamingSnapshotManager {
    /// Create a new snapshot manager.
    pub fn new() -> Self {
        Self::default()
    }

    /// Take a snapshot of the current naming context.
    pub fn take_snapshot(
        &mut self,
        context: &EnhancedNamingContext,
        operation: Option<String>,
    ) -> u64 {
        // Truncate any redo history.
        self.snapshots.truncate(self.current_index + 1);

        let id = self.next_id;
        self.next_id += 1;

        let snapshot = NamingContextSnapshot::from_context(context, id, operation);
        self.snapshots.push(snapshot);
        self.current_index = self.snapshots.len() - 1;

        id
    }

    /// Check if undo is available.
    pub fn can_undo(&self) -> bool {
        self.current_index > 0
    }

    /// Check if redo is available.
    pub fn can_redo(&self) -> bool {
        self.current_index + 1 < self.snapshots.len()
    }

    /// Undo to the previous snapshot.
    pub fn undo(&mut self, context: &mut EnhancedNamingContext) -> Option<&NamingContextSnapshot> {
        if !self.can_undo() {
            return None;
        }
        self.current_index -= 1;
        self.snapshots[self.current_index].restore_to(context);
        Some(&self.snapshots[self.current_index])
    }

    /// Redo to the next snapshot.
    pub fn redo(&mut self, context: &mut EnhancedNamingContext) -> Option<&NamingContextSnapshot> {
        if !self.can_redo() {
            return None;
        }
        self.current_index += 1;
        self.snapshots[self.current_index].restore_to(context);
        Some(&self.snapshots[self.current_index])
    }

    /// Get the current snapshot.
    pub fn current(&self) -> Option<&NamingContextSnapshot> {
        self.snapshots.get(self.current_index)
    }

    /// Get the number of snapshots.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    /// Check if there are no snapshots.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    /// Clear all snapshots.
    pub fn clear(&mut self) {
        self.snapshots.clear();
        self.current_index = 0;
    }

    /// Get the undo history (snapshots before current).
    pub fn undo_history(&self) -> &[NamingContextSnapshot] {
        &self.snapshots[..self.current_index]
    }

    /// Get the redo history (snapshots after current).
    pub fn redo_history(&self) -> &[NamingContextSnapshot] {
        &self.snapshots[self.current_index + 1..]
    }
}
