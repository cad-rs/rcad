//! Assembly / instancing scene tree.
//!
//! # Design
//!
//! - **`Arc<BRep>` sharing**: one body can be referenced from many [`AssemblyNode`]s without copying mesh data.
//! - **`DAffine3` transforms**: each node stores its transform relative to its parent; queries accumulate along the path.
//! - **Two expansion modes**:
//!   - [`Assembly::flatten`] — `(Arc<BRep>, world transform)` per leaf (lazy instancing).
//!   - [`Assembly::to_brep`] — fuse into one `BRep` via `BRep::transformed` + `append_brep`.

use std::collections::BTreeMap;
use std::sync::Arc;

use glam::{DAffine3, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::topods;
use serde::{Deserialize, Serialize};

use crate::append_brep;

/// Semantic metadata carried by assembly/document nodes.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AssemblyMetadata {
    /// Optional display label.
    pub display_name: Option<String>,
    /// Optional layer tag.
    pub layer: Option<String>,
    /// Optional material tag.
    pub material: Option<String>,
    /// Free-form key-value attributes.
    pub attributes: BTreeMap<String, String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Core types
// ─────────────────────────────────────────────────────────────────────────────

/// One node in the assembly tree.
///
/// Either a leaf holding a concrete [`BRep`] or a sub-assembly of child nodes.
/// `transform` is relative to the parent (defaults to identity).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyNode {
    /// Stable id assigned by the owning [`Assembly`].
    pub id: u64,
    /// Human-readable label.
    pub name: String,
    /// Affine transform relative to the parent (translation / rotation / scale).
    pub transform: DAffine3,
    /// Payload: shared geometry or nested nodes.
    pub content: NodeContent,
    /// Semantic metadata for this node.
    #[serde(default)]
    pub metadata: AssemblyMetadata,
}

/// Payload stored inside an [`AssemblyNode`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum NodeContent {
    /// Leaf referencing shared `BRep` data.
    Leaf(Arc<BRep>),
    /// Branch node containing an ordered list of children.
    Assembly(Vec<AssemblyNode>),
}

/// Root of an assembly document: zero or more top-level [`AssemblyNode`]s.
///
/// # Example
///
/// ```
/// # use std::sync::Arc;
/// # use glam::{DAffine3, DVec3};
/// # use rcad_kernel::BRep;
/// use rcad_scene::assembly::Assembly;
///
/// let box_brep = Arc::new(BRep::new());
/// let mut asm = Assembly::new("my_assembly");
/// asm.add_part("part_a", Arc::clone(&box_brep));
/// asm.add_part_with_transform(
///     "part_b",
///     Arc::clone(&box_brep),
///     DAffine3::from_translation(DVec3::new(5.0, 0.0, 0.0)),
/// );
/// let flat = asm.flatten();
/// assert_eq!(flat.len(), 2);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Assembly {
    /// Document / assembly name.
    pub name: String,
    /// Top-level roots (parts or nested assemblies).
    pub roots: Vec<AssemblyNode>,
    /// Metadata for the whole assembly document.
    #[serde(default)]
    pub metadata: AssemblyMetadata,
    /// Monotonic id generator for new nodes.
    next_id: u64,
}

// ─────────────────────────────────────────────────────────────────────────────
// `AssemblyNode` impl
// ─────────────────────────────────────────────────────────────────────────────

impl AssemblyNode {
    /// Leaf with identity transform.
    pub fn new_leaf(id: u64, name: impl Into<String>, brep: Arc<BRep>) -> Self {
        Self {
            id,
            name: name.into(),
            transform: DAffine3::IDENTITY,
            content: NodeContent::Leaf(brep),
            metadata: AssemblyMetadata::default(),
        }
    }

    /// Leaf with an explicit parent-relative transform.
    pub fn new_leaf_with_transform(
        id: u64,
        name: impl Into<String>,
        brep: Arc<BRep>,
        transform: DAffine3,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            transform,
            content: NodeContent::Leaf(brep),
            metadata: AssemblyMetadata::default(),
        }
    }

    /// Empty branch node (identity transform, no children yet).
    pub fn new_assembly(id: u64, name: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            transform: DAffine3::IDENTITY,
            content: NodeContent::Assembly(Vec::new()),
            metadata: AssemblyMetadata::default(),
        }
    }

    /// DFS flatten: append every leaf under this node using `parent_xform * self.transform`.
    fn flatten_into(&self, parent_xform: DAffine3, out: &mut Vec<(Arc<BRep>, DAffine3)>) {
        let world = parent_xform * self.transform;
        match &self.content {
            NodeContent::Leaf(brep) => {
                out.push((Arc::clone(brep), world));
            }
            NodeContent::Assembly(children) => {
                for child in children {
                    child.flatten_into(world, out);
                }
            }
        }
    }

    /// Recursively merge leaves into `dst`, applying accumulated transforms.
    fn merge_into(&self, parent_xform: DAffine3, dst: &mut BRep) {
        let world = parent_xform * self.transform;
        match &self.content {
            NodeContent::Leaf(brep) => {
                let mut materialized = (**brep).clone();
                materialized.apply_transform(world);
                append_brep(dst, materialized);
            }
            NodeContent::Assembly(children) => {
                for child in children {
                    child.merge_into(world, dst);
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// `Assembly` impl
// ─────────────────────────────────────────────────────────────────────────────

impl Assembly {
    /// Create an empty assembly shell.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            roots: Vec::new(),
            metadata: AssemblyMetadata::default(),
            next_id: 1,
        }
    }

    /// Attach or overwrite a document-level attribute.
    pub fn set_attribute(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.metadata.attributes.insert(key.into(), value.into());
    }

    /// Set assembly-level layer tag.
    pub fn set_layer(&mut self, layer: impl Into<String>) {
        self.metadata.layer = Some(layer.into());
    }

    /// Set assembly-level material tag.
    pub fn set_material(&mut self, material: impl Into<String>) {
        self.metadata.material = Some(material.into());
    }

    // --- internal id allocation ---

    fn alloc_id(&mut self) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        id
    }

    // --- add nodes ---

    /// Append a root leaf (identity transform); returns the new id.
    pub fn add_part(&mut self, name: impl Into<String>, brep: Arc<BRep>) -> u64 {
        let id = self.alloc_id();
        self.roots.push(AssemblyNode::new_leaf(id, name, brep));
        id
    }

    /// Append a root leaf with transform; returns the new id.
    pub fn add_part_with_transform(
        &mut self,
        name: impl Into<String>,
        brep: Arc<BRep>,
        transform: DAffine3,
    ) -> u64 {
        let id = self.alloc_id();
        self.roots.push(AssemblyNode::new_leaf_with_transform(
            id, name, brep, transform,
        ));
        id
    }

    /// Insert a fully built root node (leaf or branch); returns its id.
    pub fn add_node(&mut self, mut node: AssemblyNode) -> u64 {
        // `id == 0` means "please assign"
        if node.id == 0 {
            node.id = self.alloc_id();
        } else {
            // Keep the allocator ahead of user-supplied ids
            if node.id >= self.next_id {
                self.next_id = node.id + 1;
            }
        }
        let id = node.id;
        self.roots.push(node);
        id
    }

    // --- queries ---

    /// Flatten into `(Arc<BRep>, world transform)` pairs.
    ///
    /// One entry per leaf instance; transforms include the full parent chain for rendering
    /// without copying `BRep` vertices.
    pub fn flatten(&self) -> Vec<(Arc<BRep>, DAffine3)> {
        let mut out = Vec::new();
        for root in &self.roots {
            root.flatten_into(DAffine3::IDENTITY, &mut out);
        }
        out
    }

    /// Fuse the hierarchy into one [`BRep`].
    ///
    /// Each leaf is copied via [`BRep::transformed`], then concatenated with [`append_brep`].
    /// Use this when downstream algorithms require a single solid (Booleans, monolithic STEP, …).
    pub fn to_brep(&self) -> BRep {
        let mut result = BRep::new();
        for root in &self.roots {
            root.merge_into(DAffine3::IDENTITY, &mut result);
        }
        result
    }

    /// Total leaf count (number of instanced bodies).
    pub fn instance_count(&self) -> usize {
        fn count(node: &AssemblyNode) -> usize {
            match &node.content {
                NodeContent::Leaf(_) => 1,
                NodeContent::Assembly(children) => children.iter().map(count).sum(),
            }
        }
        self.roots.iter().map(count).sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers for `rcad_step::AssemblyComponent`
// ─────────────────────────────────────────────────────────────────────────────

/// Build an [`Assembly`] from flat `(name, brep, transform)` tuples.
///
/// Convenient when adapting `rcad_step::read_assembly` output into the scene tree.
pub fn assembly_from_parts(
    name: impl Into<String>,
    parts: impl IntoIterator<Item = (String, BRep, DAffine3)>,
) -> Assembly {
    let mut asm = Assembly::new(name);
    for (part_name, brep, transform) in parts {
        asm.add_part_with_transform(part_name, Arc::new(brep), transform);
    }
    asm
}

// ─────────────────────────────────────────────────────────────────────────────
// `DVec3` convenience helpers
// ─────────────────────────────────────────────────────────────────────────────

impl Assembly {
    /// Append a translated root leaf (rotation/scale stay identity).
    pub fn add_part_at(
        &mut self,
        name: impl Into<String>,
        brep: Arc<BRep>,
        translation: DVec3,
    ) -> u64 {
        self.add_part_with_transform(name, brep, DAffine3::from_translation(translation))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembly_metadata_setters_work() {
        let mut asm = Assembly::new("m");
        asm.set_layer("L1");
        asm.set_material("Aluminum");
        asm.set_attribute("owner", "team-a");

        assert_eq!(asm.metadata.layer.as_deref(), Some("L1"));
        assert_eq!(asm.metadata.material.as_deref(), Some("Aluminum"));
        assert_eq!(
            asm.metadata.attributes.get("owner").map(String::as_str),
            Some("team-a")
        );
    }

    #[test]
    fn node_metadata_defaults_empty() {
        let node = AssemblyNode::new_assembly(1, "sub");
        assert!(node.metadata.layer.is_none());
        assert!(node.metadata.attributes.is_empty());
    }
}
