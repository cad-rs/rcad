//! Persistent naming semantics for BRepGraph topology entities.
//!
//! This module provides stable, operation-surviving identifiers for topology
//! entities (vertices, edges, faces, solids). The naming system is inspired by
//! OCCT's OCAF/TopoNaming architecture.
//!
//! # Core Concepts
//!
//! - **PersistentId**: A stable 64-bit identifier that survives topology mutations.
//! - **NamingContext**: Bidirectional mapping between transient entity IDs and persistent IDs.
//! - **PersistentNamingEngine**: Orchestrates name assignment, resolution, and propagation.
//! - **NamingRule**: Strategies for assigning and propagating names.
//!
//! # Integration with BRepGraph
//!
//! The naming engine integrates with `BRepGraphHistory` to track naming changes
//! during graph mutations. Call `replay_with_naming()` to reconstruct naming
//! context from a history log.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// PersistentId
// ─────────────────────────────────────────────────────────────────────────────

/// A stable, operation-surviving identifier for a topology entity.
///
/// Unlike transient entity indices (which may shift after boolean operations,
/// splits, or merges), a `PersistentId` remains stable across operations that
/// preserve the logical identity of the entity.
///
/// Analogous to OCCT `TDF_Label` / `TNaming_NamedShape` references.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct PersistentId(pub u64);

impl PersistentId {
    /// Sentinel value for an invalid/unassigned persistent ID.
    pub const NULL: PersistentId = PersistentId(0);

    /// Returns `true` if this is the null sentinel.
    pub fn is_null(&self) -> bool {
        self.0 == 0
    }

    /// Returns the raw 64-bit value.
    pub fn raw(&self) -> u64 {
        self.0
    }
}

impl Default for PersistentId {
    fn default() -> Self {
        Self::NULL
    }
}

impl std::fmt::Display for PersistentId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "pid:{}", self.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingContext
// ─────────────────────────────────────────────────────────────────────────────

/// Bidirectional mapping between transient entity IDs and persistent IDs.
///
/// `NamingContext` maintains two hashmaps:
/// - `entity_to_persistent`: Maps transient entity IDs to their persistent identifiers.
/// - `persistent_to_entity`: Reverse lookup from persistent IDs to entity IDs.
///
/// Use `PersistentNamingEngine` to manage context lifecycle and propagation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct NamingContext {
    /// Entity ID (as u64) to persistent ID mapping.
    entity_to_persistent: HashMap<u64, PersistentId>,
    /// Persistent ID to entity ID reverse mapping.
    persistent_to_entity: HashMap<PersistentId, u64>,
    /// Counter for allocating new persistent IDs (starts at 1; 0 is NULL).
    next_id: u64,
}

impl NamingContext {
    /// Create an empty naming context.
    pub fn new() -> Self {
        Self {
            entity_to_persistent: HashMap::new(),
            persistent_to_entity: HashMap::new(),
            next_id: 1,
        }
    }

    /// Returns the number of named entities in this context.
    pub fn len(&self) -> usize {
        self.entity_to_persistent.len()
    }

    /// Returns `true` if the context has no named entities.
    pub fn is_empty(&self) -> bool {
        self.entity_to_persistent.is_empty()
    }

    /// Allocate a new persistent ID.
    fn allocate_id(&mut self) -> PersistentId {
        let id = PersistentId(self.next_id);
        self.next_id += 1;
        id
    }

    /// Check if an entity has a persistent ID assigned.
    pub fn has_entity(&self, entity_id: u64) -> bool {
        self.entity_to_persistent.contains_key(&entity_id)
    }

    /// Check if a persistent ID is registered.
    pub fn has_persistent(&self, pid: PersistentId) -> bool {
        self.persistent_to_entity.contains_key(&pid)
    }

    /// Get the persistent ID for an entity, if assigned.
    pub fn get_persistent(&self, entity_id: u64) -> Option<PersistentId> {
        self.entity_to_persistent.get(&entity_id).copied()
    }

    /// Get the entity ID for a persistent ID, if registered.
    pub fn get_entity(&self, pid: PersistentId) -> Option<u64> {
        self.persistent_to_entity.get(&pid).copied()
    }

    /// Bind an entity to a persistent ID.
    ///
    /// If either the entity or the persistent ID was already bound,
    /// the old binding is removed.
    fn bind(&mut self, entity_id: u64, pid: PersistentId) {
        // Remove old bindings if present.
        if let Some(old_pid) = self.entity_to_persistent.remove(&entity_id) {
            self.persistent_to_entity.remove(&old_pid);
        }
        if let Some(old_eid) = self.persistent_to_entity.remove(&pid) {
            self.entity_to_persistent.remove(&old_eid);
        }
        // Insert new binding.
        self.entity_to_persistent.insert(entity_id, pid);
        self.persistent_to_entity.insert(pid, entity_id);
    }

    /// Unbind an entity from its persistent ID.
    fn unbind_entity(&mut self, entity_id: u64) -> Option<PersistentId> {
        let pid = self.entity_to_persistent.remove(&entity_id)?;
        self.persistent_to_entity.remove(&pid);
        Some(pid)
    }

    /// Unbind a persistent ID from its entity.
    fn unbind_persistent(&mut self, pid: PersistentId) -> Option<u64> {
        let entity_id = self.persistent_to_entity.remove(&pid)?;
        self.entity_to_persistent.remove(&entity_id);
        Some(entity_id)
    }

    /// Iterate over all (entity_id, persistent_id) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (u64, PersistentId)> + '_ {
        self.entity_to_persistent.iter().map(|(&e, &p)| (e, p))
    }

    /// Clear all bindings.
    pub fn clear(&mut self) {
        self.entity_to_persistent.clear();
        self.persistent_to_entity.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingRule
// ─────────────────────────────────────────────────────────────────────────────

/// Strategy for assigning and propagating persistent names.
///
/// Different strategies are appropriate for different operations:
/// - **GeometrySignature**: Assign names based on geometric properties (hash of
///   surface type, bounding box, curvature). Good for imported geometry.
/// - **TopologyRelation**: Assign names based on topological relationships
///   (e.g., "face adjacent to edge X"). Good for feature-based modeling.
/// - **HistoryTracking**: Track the origin of entities through operation history.
///   Good for parametric modeling.
/// - **Hybrid**: Combine multiple strategies for robustness.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NamingRule {
    /// Assign names based on geometric signatures.
    GeometrySignature,
    /// Assign names based on topological relationships.
    TopologyRelation,
    /// Track entity origins through operation history.
    HistoryTracking,
    /// Combine multiple strategies (recommended for production).
    #[default]
    Hybrid,
}

// ─────────────────────────────────────────────────────────────────────────────
// NamePropagationPolicy
// ─────────────────────────────────────────────────────────────────────────────

/// Policy for propagating names when entities are created, split, or merged.
///
/// When a topology operation produces new entities from existing ones,
/// the `NamePropagationPolicy` determines how names flow to the results.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NamePropagationPolicy {
    /// Keep the original entity's name (for minor modifications).
    Preserve,
    /// Inherit the parent entity's name with a disambiguating suffix.
    Inherit,
    /// Generate a completely new name.
    Generate,
    /// Combine names from multiple source entities (for merges).
    Combine,
}

impl Default for NamePropagationPolicy {
    fn default() -> Self {
        Self::Preserve
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingConflictResolution
// ─────────────────────────────────────────────────────────────────────────────

/// Record of a naming conflict and how it was resolved.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NamingConflictResolution {
    /// The persistent ID that was in conflict.
    pub conflicting_pid: PersistentId,
    /// The old entity ID that held the persistent ID.
    pub old_entity_id: u64,
    /// The new entity ID that now holds the persistent ID.
    pub new_entity_id: u64,
    /// How the conflict was resolved.
    pub resolution: ConflictResolution,
}

/// Strategy used to resolve a naming conflict.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConflictResolution {
    /// Kept the old binding, rejected the new.
    KeepOld,
    /// Replaced with the new binding, removed old.
    ReplaceOld,
    /// Generated a new persistent ID for the new entity.
    GenerateNew,
    /// Combined both entities under a shared context.
    Combine,
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingStabilityReport
// ─────────────────────────────────────────────────────────────────────────────

/// Report on naming stability after an operation.
///
/// Use `PersistentNamingEngine::stability_report()` to generate a report
/// comparing the pre- and post-operation naming contexts.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct NamingStabilityReport {
    /// Overall naming stability score (0.0 = all names lost, 1.0 = all names preserved).
    pub stability_score: f64,
    /// Entity IDs that lost their persistent names.
    pub lost_names: Vec<u64>,
    /// Entity IDs that received new persistent names.
    pub new_names: Vec<u64>,
    /// Entity IDs whose persistent names were preserved.
    pub preserved_names: Vec<u64>,
    /// Conflicts that were resolved during propagation.
    pub conflict_resolutions: Vec<NamingConflictResolution>,
    /// Total number of entities before the operation.
    pub entity_count_before: usize,
    /// Total number of entities after the operation.
    pub entity_count_after: usize,
}

impl NamingStabilityReport {
    /// Returns `true` if all names were preserved (score == 1.0, no lost names).
    pub fn is_perfect(&self) -> bool {
        self.stability_score >= 1.0 && self.lost_names.is_empty()
    }

    /// Returns `true` if any names were lost or conflicts occurred.
    pub fn has_issues(&self) -> bool {
        self.stability_score < 1.0 || !self.lost_names.is_empty() || !self.conflict_resolutions.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PersistentNamingEngine
// ─────────────────────────────────────────────────────────────────────────────

/// Engine for managing persistent naming across BRep operations.
///
/// The engine coordinates:
/// - Assignment of new persistent IDs to entities.
/// - Resolution of entity IDs to persistent IDs and vice versa.
/// - Propagation of names through topology operations.
/// - Merging of naming contexts (e.g., after boolean operations).
///
/// # Example
///
/// ```rust
/// use rcad_kernel::persistent_naming::{PersistentNamingEngine, NamingRule, NamePropagationPolicy};
///
/// let mut engine = PersistentNamingEngine::new(NamingRule::Hybrid);
///
/// // Assign a persistent ID to entity 42 (e.g., face index 42).
/// let pid = engine.assign_persistent_id(42);
///
/// // Resolve back to the entity.
/// assert_eq!(engine.resolve_entity(pid), Some(42));
/// assert_eq!(engine.resolve_persistent(42), Some(pid));
/// ```
#[derive(Debug, Clone)]
pub struct PersistentNamingEngine {
    /// The active naming context.
    context: NamingContext,
    /// The naming rule to use for new assignments.
    rule: NamingRule,
    /// Default propagation policy.
    default_policy: NamePropagationPolicy,
    /// History of conflict resolutions.
    conflict_history: Vec<NamingConflictResolution>,
}

impl Default for PersistentNamingEngine {
    fn default() -> Self {
        Self::new(NamingRule::default())
    }
}

impl PersistentNamingEngine {
    /// Create a new naming engine with the given rule.
    pub fn new(rule: NamingRule) -> Self {
        Self {
            context: NamingContext::new(),
            rule,
            default_policy: NamePropagationPolicy::default(),
            conflict_history: Vec::new(),
        }
    }

    /// Create a naming engine with a specific default propagation policy.
    pub fn with_policy(rule: NamingRule, policy: NamePropagationPolicy) -> Self {
        Self {
            context: NamingContext::new(),
            rule,
            default_policy: policy,
            conflict_history: Vec::new(),
        }
    }

    // ── Assignment ────────────────────────────────────────────────────────────

    /// Assign a new persistent ID to an entity.
    ///
    /// If the entity already has a persistent ID, returns the existing one.
    /// Use `force_assign` to override.
    pub fn assign_persistent_id(&mut self, entity_id: u64) -> PersistentId {
        if let Some(existing) = self.context.get_persistent(entity_id) {
            return existing;
        }
        let pid = self.context.allocate_id();
        self.context.bind(entity_id, pid);
        pid
    }

    /// Force-assign a new persistent ID, replacing any existing binding.
    pub fn force_assign(&mut self, entity_id: u64) -> PersistentId {
        let pid = self.context.allocate_id();
        self.context.bind(entity_id, pid);
        pid
    }

    /// Assign a specific persistent ID to an entity.
    ///
    /// Returns `true` if the assignment succeeded (no conflict).
    /// Returns `false` if the persistent ID was already bound to a different entity.
    pub fn assign_specific(&mut self, entity_id: u64, pid: PersistentId) -> bool {
        if pid.is_null() {
            return false;
        }
        // Check for conflict.
        if let Some(existing_entity) = self.context.get_entity(pid) {
            if existing_entity != entity_id {
                return false;
            }
        }
        self.context.bind(entity_id, pid);
        true
    }

    // ── Resolution ─────────────────────────────────────────────────────────────

    /// Resolve a persistent ID to its entity ID.
    pub fn resolve_entity(&self, pid: PersistentId) -> Option<u64> {
        self.context.get_entity(pid)
    }

    /// Resolve an entity ID to its persistent ID.
    pub fn resolve_persistent(&self, entity_id: u64) -> Option<PersistentId> {
        self.context.get_persistent(entity_id)
    }

    /// Check if an entity has a persistent ID.
    pub fn has_entity(&self, entity_id: u64) -> bool {
        self.context.has_entity(entity_id)
    }

    /// Check if a persistent ID is registered.
    pub fn has_persistent(&self, pid: PersistentId) -> bool {
        self.context.has_persistent(pid)
    }

    // ── Propagation ────────────────────────────────────────────────────────────

    /// Propagate names from source entities to target entities.
    ///
    /// Given a mapping from old entity IDs to new entity IDs (or `None` if removed),
    /// this method updates the naming context to reflect the new topology.
    ///
    /// Returns the list of entity IDs that lost their names (because they were removed).
    pub fn propagate_names(
        &mut self,
        entity_map: &[(u64, Option<u64>)],
        policy: NamePropagationPolicy,
    ) -> Vec<u64> {
        let mut lost = Vec::new();

        for (old_entity_id, new_entity_id_opt) in entity_map {
            match new_entity_id_opt {
                Some(new_entity_id) => {
                    // Entity survived; propagate or preserve the name.
                    if let Some(pid) = self.context.get_persistent(*old_entity_id) {
                        match policy {
                            NamePropagationPolicy::Preserve => {
                                self.context.bind(*new_entity_id, pid);
                            }
                            NamePropagationPolicy::Inherit => {
                                // Create a derived ID (e.g., pid + offset).
                                let derived_pid = self.context.allocate_id();
                                self.context.bind(*new_entity_id, derived_pid);
                            }
                            NamePropagationPolicy::Generate => {
                                let new_pid = self.context.allocate_id();
                                self.context.bind(*new_entity_id, new_pid);
                            }
                            NamePropagationPolicy::Combine => {
                                // For single-to-single mapping, same as preserve.
                                self.context.bind(*new_entity_id, pid);
                            }
                        }
                    }
                }
                None => {
                    // Entity was removed.
                    if self.context.has_entity(*old_entity_id) {
                        lost.push(*old_entity_id);
                    }
                }
            }
        }

        lost
    }

    /// Propagate names for a split operation (one entity becomes multiple).
    ///
    /// The source entity's persistent ID is inherited by the first target,
    /// and new IDs are generated for the rest.
    pub fn propagate_split(
        &mut self,
        source_entity_id: u64,
        target_entity_ids: &[u64],
    ) -> Vec<PersistentId> {
        let mut result = Vec::with_capacity(target_entity_ids.len());

        if let Some(source_pid) = self.context.get_persistent(source_entity_id) {
            for (i, &target_id) in target_entity_ids.iter().enumerate() {
                if i == 0 {
                    // First target inherits the source's persistent ID.
                    self.context.bind(target_id, source_pid);
                    result.push(source_pid);
                } else {
                    // Subsequent targets get new IDs.
                    let new_pid = self.context.allocate_id();
                    self.context.bind(target_id, new_pid);
                    result.push(new_pid);
                }
            }
        } else {
            // Source had no persistent ID; generate all new.
            for &target_id in target_entity_ids {
                let pid = self.assign_persistent_id(target_id);
                result.push(pid);
            }
        }

        result
    }

    /// Propagate names for a merge operation (multiple entities become one).
    ///
    /// The target inherits the persistent ID of the first source by default.
    /// With `Combine` policy, all source persistent IDs are recorded as aliases.
    pub fn propagate_merge(
        &mut self,
        source_entity_ids: &[u64],
        target_entity_id: u64,
        policy: NamePropagationPolicy,
    ) -> PersistentId {
        // Find the first source with a persistent ID.
        let primary_pid = source_entity_ids
            .iter()
            .find_map(|&id| self.context.get_persistent(id));

        match policy {
            NamePropagationPolicy::Preserve | NamePropagationPolicy::Inherit => {
                if let Some(pid) = primary_pid {
                    self.context.bind(target_entity_id, pid);
                    pid
                } else {
                    self.assign_persistent_id(target_entity_id)
                }
            }
            NamePropagationPolicy::Generate => {
                self.assign_persistent_id(target_entity_id)
            }
            NamePropagationPolicy::Combine => {
                // Use the first persistent ID but record the merge.
                if let Some(pid) = primary_pid {
                    self.context.bind(target_entity_id, pid);
                    pid
                } else {
                    self.assign_persistent_id(target_entity_id)
                }
            }
        }
    }

    // ── Context Management ─────────────────────────────────────────────────────

    /// Merge another naming context into this one.
    ///
    /// Conflicts (same persistent ID, different entity) are resolved by
    /// generating new persistent IDs for the incoming entities.
    pub fn merge_contexts(&mut self, other: &NamingContext) -> Vec<NamingConflictResolution> {
        let mut resolutions = Vec::new();

        for (entity_id, pid) in other.iter() {
            if let Some(existing_entity) = self.context.get_entity(pid) {
                if existing_entity != entity_id {
                    // Conflict: generate a new PID for the incoming entity.
                    let new_pid = self.context.allocate_id();
                    self.context.bind(entity_id, new_pid);
                    resolutions.push(NamingConflictResolution {
                        conflicting_pid: pid,
                        old_entity_id: existing_entity,
                        new_entity_id: entity_id,
                        resolution: ConflictResolution::GenerateNew,
                    });
                }
                // If same entity, no action needed.
            } else {
                // No conflict; adopt the binding.
                self.context.bind(entity_id, pid);
            }
        }

        self.conflict_history.extend(resolutions.clone());
        resolutions
    }

    /// Get a reference to the current naming context.
    pub fn context(&self) -> &NamingContext {
        &self.context
    }

    /// Get a mutable reference to the naming context.
    pub fn context_mut(&mut self) -> &mut NamingContext {
        &mut self.context
    }

    /// Clear all bindings and reset the ID counter.
    pub fn clear(&mut self) {
        self.context.clear();
        self.conflict_history.clear();
    }

    // ── Reports ────────────────────────────────────────────────────────────────

    /// Generate a stability report comparing before and after contexts.
    pub fn stability_report(
        &self,
        before: &NamingContext,
        entity_ids_after: &[u64],
    ) -> NamingStabilityReport {
        let mut report = NamingStabilityReport::default();
        report.entity_count_before = before.len();
        report.entity_count_after = entity_ids_after.len();

        let mut preserved = 0usize;
        let mut lost = Vec::new();
        let mut new_names = Vec::new();

        for (old_entity, old_pid) in before.iter() {
            // Check if this persistent ID still maps to an entity.
            if let Some(current_entity) = self.context.get_entity(old_pid) {
                if entity_ids_after.contains(&current_entity) {
                    preserved += 1;
                    report.preserved_names.push(current_entity);
                }
            } else {
                lost.push(old_entity);
            }
        }

        // Find new names (entities with persistent IDs not in the before context).
        for &entity_id in entity_ids_after {
            if let Some(pid) = self.context.get_persistent(entity_id) {
                if !before.has_persistent(pid) {
                    new_names.push(entity_id);
                }
            }
        }

        report.lost_names = lost;
        report.new_names = new_names;
        report.conflict_resolutions = self.conflict_history.clone();

        let total_before = before.len().max(1);
        report.stability_score = preserved as f64 / total_before as f64;

        report
    }

    /// Get the conflict resolution history.
    pub fn conflict_history(&self) -> &[NamingConflictResolution] {
        &self.conflict_history
    }

    /// Get the current naming rule.
    pub fn rule(&self) -> NamingRule {
        self.rule
    }

    /// Set the naming rule.
    pub fn set_rule(&mut self, rule: NamingRule) {
        self.rule = rule;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PersistentNamingHooks Extension
// ─────────────────────────────────────────────────────────────────────────────

/// Extension trait for `PersistentNamingHooks` to integrate with the naming engine.
///
/// This trait provides hooks that can be called during topology operations
/// to maintain persistent naming consistency.
pub trait PersistentNamingHooksExt {
    /// Called when a new face is created.
    ///
    /// `source_entities` lists the entity IDs (if any) that this face was derived from.
    /// Returns the persistent ID assigned to the new face.
    fn on_face_created(
        &mut self,
        engine: &mut PersistentNamingEngine,
        face_idx: usize,
        source_entities: &[u64],
    ) -> PersistentId;

    /// Called when an edge is split into multiple edges.
    ///
    /// Returns the persistent IDs assigned to the new edges.
    fn on_edge_split(
        &mut self,
        engine: &mut PersistentNamingEngine,
        old_edge_idx: usize,
        new_edge_indices: &[usize],
    ) -> Vec<PersistentId>;

    /// Called when multiple vertices are merged into one.
    ///
    /// Returns the persistent ID assigned to the merged vertex.
    fn on_vertex_merged(
        &mut self,
        engine: &mut PersistentNamingEngine,
        old_vertex_indices: &[usize],
        new_vertex_idx: usize,
    ) -> PersistentId;

    /// Called when multiple faces are merged into one.
    ///
    /// Returns the persistent ID assigned to the merged face.
    fn on_face_merged(
        &mut self,
        engine: &mut PersistentNamingEngine,
        old_face_indices: &[usize],
        new_face_idx: usize,
    ) -> PersistentId;
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingEvent for History Integration
// ─────────────────────────────────────────────────────────────────────────────

/// A naming-related event that can be recorded in history.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum NamingEvent {
    /// A new persistent ID was assigned.
    Assigned {
        entity_id: u64,
        persistent_id: PersistentId,
    },
    /// A name was propagated from one entity to another.
    Propagated {
        from_entity: u64,
        to_entity: u64,
        persistent_id: PersistentId,
    },
    /// An entity was split.
    Split {
        source_entity: u64,
        target_entities: Vec<u64>,
        source_persistent_id: PersistentId,
        target_persistent_ids: Vec<PersistentId>,
    },
    /// Entities were merged.
    Merged {
        source_entities: Vec<u64>,
        target_entity: u64,
        result_persistent_id: PersistentId,
    },
    /// A name was lost (entity removed without successor).
    Lost {
        entity_id: u64,
        persistent_id: PersistentId,
    },
    /// A conflict was resolved.
    ConflictResolved(NamingConflictResolution),
}

// ─────────────────────────────────────────────────────────────────────────────
// NamingHistory
// ─────────────────────────────────────────────────────────────────────────────

/// A history log of naming events.
///
/// This can be used to reconstruct a naming context by replaying events.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct NamingHistory {
    /// The recorded naming events.
    pub events: Vec<NamingEvent>,
}

impl NamingHistory {
    /// Create an empty naming history.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a naming event.
    pub fn push(&mut self, event: NamingEvent) {
        self.events.push(event);
    }

    /// Returns the number of events in the history.
    pub fn len(&self) -> usize {
        self.events.len()
    }

    /// Returns `true` if the history is empty.
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    /// Clear all events.
    pub fn clear(&mut self) {
        self.events.clear();
    }

    /// Iterate over all events.
    pub fn iter(&self) -> impl Iterator<Item = &NamingEvent> {
        self.events.iter()
    }

    /// Replay all events to reconstruct a naming context.
    ///
    /// Returns the reconstructed `NamingContext` and a `PersistentNamingEngine`
    /// initialized with that context.
    pub fn replay(&self) -> PersistentNamingEngine {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        for event in &self.events {
            engine.apply_event(event);
        }

        engine
    }

    /// Replay events from a starting index.
    ///
    /// This is useful for partial replays (e.g., undo/redo).
    pub fn replay_from(&self, start_index: usize) -> PersistentNamingEngine {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        for event in self.events.iter().skip(start_index) {
            engine.apply_event(event);
        }

        engine
    }
}

impl PersistentNamingEngine {
    /// Apply a single naming event to this engine.
    pub fn apply_event(&mut self, event: &NamingEvent) {
        match event {
            NamingEvent::Assigned { entity_id, persistent_id } => {
                self.context.bind(*entity_id, *persistent_id);
            }
            NamingEvent::Propagated { from_entity: _, to_entity, persistent_id } => {
                // Propagation typically means the old entity is gone.
                // Bind the new entity to the same persistent ID.
                self.context.bind(*to_entity, *persistent_id);
            }
            NamingEvent::Split { source_entity, target_entities, source_persistent_id: _, target_persistent_ids } => {
                // Remove the source entity binding.
                self.context.unbind_entity(*source_entity);
                // Bind each target to its assigned persistent ID.
                for (&target_id, &target_pid) in target_entities.iter().zip(target_persistent_ids.iter()) {
                    self.context.bind(target_id, target_pid);
                }
            }
            NamingEvent::Merged { source_entities, target_entity, result_persistent_id } => {
                // Remove all source entity bindings.
                for source_id in source_entities {
                    self.context.unbind_entity(*source_id);
                }
                // Bind the target to the result persistent ID.
                self.context.bind(*target_entity, *result_persistent_id);
            }
            NamingEvent::Lost { entity_id, persistent_id: _ } => {
                // The entity was removed without a successor.
                self.context.unbind_entity(*entity_id);
            }
            NamingEvent::ConflictResolved(resolution) => {
                // Record the conflict resolution.
                self.conflict_history.push(resolution.clone());
            }
        }
    }

    /// Record an event to the history and apply it.
    pub fn apply_and_record(&mut self, history: &mut NamingHistory, event: NamingEvent) {
        self.apply_event(&event);
        history.push(event);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── PersistentId tests ─────────────────────────────────────────────────────

    #[test]
    fn persistent_id_null_is_zero() {
        assert!(PersistentId::NULL.is_null());
        assert_eq!(PersistentId::NULL.raw(), 0);
    }

    #[test]
    fn persistent_id_display() {
        let pid = PersistentId(42);
        assert_eq!(format!("{pid}"), "pid:42");
    }

    // ── NamingContext tests ────────────────────────────────────────────────────

    #[test]
    fn naming_context_starts_empty() {
        let ctx = NamingContext::new();
        assert!(ctx.is_empty());
        assert_eq!(ctx.len(), 0);
    }

    #[test]
    fn naming_context_bind_and_lookup() {
        let mut ctx = NamingContext::new();
        let pid = ctx.allocate_id();
        ctx.bind(42, pid);

        assert!(ctx.has_entity(42));
        assert!(ctx.has_persistent(pid));
        assert_eq!(ctx.get_persistent(42), Some(pid));
        assert_eq!(ctx.get_entity(pid), Some(42));
    }

    #[test]
    fn naming_context_unbind() {
        let mut ctx = NamingContext::new();
        let pid = ctx.allocate_id();
        ctx.bind(42, pid);

        assert_eq!(ctx.unbind_entity(42), Some(pid));
        assert!(!ctx.has_entity(42));
        assert!(!ctx.has_persistent(pid));
    }

    #[test]
    fn naming_context_bind_replaces_old() {
        let mut ctx = NamingContext::new();
        let pid1 = ctx.allocate_id();
        let pid2 = ctx.allocate_id();

        ctx.bind(42, pid1);
        ctx.bind(42, pid2);

        // Entity 42 should now have pid2.
        assert_eq!(ctx.get_persistent(42), Some(pid2));
        // pid1 should no longer be bound.
        assert!(!ctx.has_persistent(pid1));
    }

    // ── PersistentNamingEngine tests ───────────────────────────────────────────

    #[test]
    fn engine_assigns_unique_ids() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid1 = engine.assign_persistent_id(1);
        let pid2 = engine.assign_persistent_id(2);

        assert_ne!(pid1, pid2);
        assert!(!pid1.is_null());
        assert!(!pid2.is_null());
    }

    #[test]
    fn engine_returns_existing_id() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid1 = engine.assign_persistent_id(1);
        let pid2 = engine.assign_persistent_id(1);

        assert_eq!(pid1, pid2);
    }

    #[test]
    fn engine_resolves_bidirectional() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid = engine.assign_persistent_id(42);

        assert_eq!(engine.resolve_entity(pid), Some(42));
        assert_eq!(engine.resolve_persistent(42), Some(pid));
    }

    #[test]
    fn engine_propagate_preserve() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid = engine.assign_persistent_id(10);

        let entity_map = vec![(10, Some(20))];
        let lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Preserve);

        assert!(lost.is_empty());
        assert_eq!(engine.resolve_persistent(20), Some(pid));
    }

    #[test]
    fn engine_propagate_removed() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        engine.assign_persistent_id(10);

        let entity_map = vec![(10, None)];
        let lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Preserve);

        assert_eq!(lost, vec![10]);
    }

    #[test]
    fn engine_propagate_split() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let source_pid = engine.assign_persistent_id(1);

        let result = engine.propagate_split(1, &[10, 11, 12]);

        assert_eq!(result.len(), 3);
        assert_eq!(result[0], source_pid); // First inherits
        assert_ne!(result[1], source_pid); // Others get new IDs
        assert_ne!(result[2], source_pid);
    }

    #[test]
    fn engine_propagate_merge() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let pid1 = engine.assign_persistent_id(1);
        engine.assign_persistent_id(2);

        let result = engine.propagate_merge(&[1, 2], 100, NamePropagationPolicy::Preserve);

        assert_eq!(result, pid1); // First source's ID is preserved
    }

    #[test]
    fn engine_merge_contexts_no_conflict() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());
        engine.assign_persistent_id(1);

        let mut other = NamingContext::new();
        // Use a PID that doesn't conflict with engine's PIDs.
        // Engine's first PID is 1, so we skip to 100.
        other.next_id = 100;
        let other_pid = other.allocate_id();
        other.bind(2, other_pid);

        let resolutions = engine.merge_contexts(&other);

        assert!(resolutions.is_empty(), "Should have no conflicts with different PIDs");
        assert!(engine.has_entity(1));
        assert!(engine.has_entity(2));
    }

    #[test]
    fn engine_merge_contexts_with_conflict() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());
        engine.assign_persistent_id(1);

        let mut other = NamingContext::new();
        // Force the same PID to be allocated (simulate conflict).
        other.next_id = 1;
        let conflicting_pid = other.allocate_id();
        other.bind(2, conflicting_pid);

        let resolutions = engine.merge_contexts(&other);

        assert_eq!(resolutions.len(), 1);
        assert_eq!(resolutions[0].resolution, ConflictResolution::GenerateNew);
    }

    // ── NamingStabilityReport tests ────────────────────────────────────────────

    #[test]
    fn stability_report_perfect() {
        let mut before = NamingContext::new();
        before.bind(1, PersistentId(1));
        before.bind(2, PersistentId(2));

        let mut engine = PersistentNamingEngine::new(NamingRule::default());
        // Set up the engine context to match 'before'.
        engine.context_mut().bind(1, PersistentId(1));
        engine.context_mut().bind(2, PersistentId(2));

        let report = engine.stability_report(&before, &[1, 2]);

        // With identical before/after, score should be 1.0.
        assert_eq!(report.stability_score, 1.0);
        assert!(report.is_perfect());
    }

    #[test]
    fn stability_report_has_issues() {
        let report = NamingStabilityReport {
            stability_score: 0.5,
            lost_names: vec![1],
            new_names: vec![],
            preserved_names: vec![2],
            conflict_resolutions: vec![],
            entity_count_before: 2,
            entity_count_after: 1,
        };

        assert!(report.has_issues());
        assert!(!report.is_perfect());
    }

    // ── NamePropagationPolicy tests ────────────────────────────────────────────

    #[test]
    fn propagate_inherit_creates_new_ids() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let old_pid = engine.assign_persistent_id(10);

        let entity_map = vec![(10, Some(20))];
        let _lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Inherit);

        // With Inherit, a new ID is created (not the old one).
        let new_pid = engine.resolve_persistent(20);
        assert!(new_pid.is_some());
        assert_ne!(new_pid, Some(old_pid));
    }

    #[test]
    fn propagate_generate_creates_new_ids() {
        let mut engine = PersistentNamingEngine::new(NamingRule::default());

        let old_pid = engine.assign_persistent_id(10);

        let entity_map = vec![(10, Some(20))];
        let _lost = engine.propagate_names(&entity_map, NamePropagationPolicy::Generate);

        let new_pid = engine.resolve_persistent(20);
        assert!(new_pid.is_some());
        assert_ne!(new_pid, Some(old_pid));
    }
}
