use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use crate::BRep;

/// A stable reference to a topological entity in a B-Rep.
///
/// Face indexing follows RCAD's flattened face order (solid/shell/face traversal).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum TopoEntityRef {
    Vertex(usize),
    Edge(usize),
    Face(usize),
    Solid(usize),
}

/// Baseline persistent naming table for topology entities.
///
/// This is a lightweight hook layer analogous to OCCT OCAF naming tables:
/// it provides stable user-level names and bidirectional resolution between
/// names and topology references.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PersistentNamingHooks {
    name_to_ref: BTreeMap<String, TopoEntityRef>,
    ref_to_name: BTreeMap<TopoEntityRef, String>,
}

impl PersistentNamingHooks {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a baseline naming table with deterministic default labels.
    ///
    /// Generated labels are:
    /// - vertices: `v0`, `v1`, ...
    /// - edges: `e0`, `e1`, ...
    /// - faces: `f0`, `f1`, ... (flattened index)
    /// - solids: `s0`, `s1`, ...
    pub fn with_default_labels_for_brep(brep: &BRep) -> Self {
        let mut out = Self::new();
        for i in 0..brep.vertices.len() {
            out.bind_unchecked(format!("v{i}"), TopoEntityRef::Vertex(i));
        }
        for i in 0..brep.edges.len() {
            out.bind_unchecked(format!("e{i}"), TopoEntityRef::Edge(i));
        }
        for i in 0..flat_face_count(brep) {
            out.bind_unchecked(format!("f{i}"), TopoEntityRef::Face(i));
        }
        for i in 0..brep.solids.len() {
            out.bind_unchecked(format!("s{i}"), TopoEntityRef::Solid(i));
        }
        out
    }

    /// Bind a user-visible `name` to an entity reference.
    ///
    /// If the name or entity is already bound, the old binding is replaced.
    pub fn bind(&mut self, name: impl Into<String>, target: TopoEntityRef) {
        self.bind_unchecked(name.into(), target);
    }

    /// Bind with topology bounds check against `brep`.
    pub fn bind_for_brep(
        &mut self,
        brep: &BRep,
        name: impl Into<String>,
        target: TopoEntityRef,
    ) -> Result<(), String> {
        if !is_valid_ref_for_brep(brep, target) {
            return Err(format!("invalid topology reference for BRep: {target:?}"));
        }
        self.bind_unchecked(name.into(), target);
        Ok(())
    }

    pub fn resolve(&self, name: &str) -> Option<TopoEntityRef> {
        self.name_to_ref.get(name).copied()
    }

    pub fn name_of(&self, target: TopoEntityRef) -> Option<&str> {
        self.ref_to_name.get(&target).map(String::as_str)
    }

    pub fn unbind_name(&mut self, name: &str) -> Option<TopoEntityRef> {
        let target = self.name_to_ref.remove(name)?;
        self.ref_to_name.remove(&target);
        Some(target)
    }

    pub fn unbind_ref(&mut self, target: TopoEntityRef) -> Option<String> {
        let name = self.ref_to_name.remove(&target)?;
        self.name_to_ref.remove(&name);
        Some(name)
    }

    pub fn rename(&mut self, old_name: &str, new_name: impl Into<String>) -> Result<(), String> {
        let Some(target) = self.resolve(old_name) else {
            return Err(format!("name '{old_name}' not found"));
        };
        let new_name = new_name.into();
        if new_name == old_name {
            return Ok(());
        }
        if let Some(existing) = self.resolve(&new_name) {
            return Err(format!(
                "name '{new_name}' is already bound to {existing:?}"
            ));
        }
        self.unbind_name(old_name);
        self.bind_unchecked(new_name, target);
        Ok(())
    }

    /// Returns all invalid bindings for the given `brep`.
    pub fn validate_against_brep(&self, brep: &BRep) -> Vec<String> {
        let mut issues = Vec::new();
        for (name, target) in &self.name_to_ref {
            if !is_valid_ref_for_brep(brep, *target) {
                issues.push(format!("name '{name}' points to out-of-range entity {target:?}"));
            }
        }
        issues
    }

    /// Remove bindings that no longer point to valid topology entities.
    pub fn retain_valid_for_brep(&mut self, brep: &BRep) {
        let invalid_names: Vec<String> = self
            .name_to_ref
            .iter()
            .filter_map(|(name, target)| {
                if is_valid_ref_for_brep(brep, *target) {
                    None
                } else {
                    Some(name.clone())
                }
            })
            .collect();
        for name in invalid_names {
            self.unbind_name(&name);
        }
    }

    pub fn len(&self) -> usize {
        self.name_to_ref.len()
    }

    pub fn is_empty(&self) -> bool {
        self.name_to_ref.is_empty()
    }

    fn bind_unchecked(&mut self, name: String, target: TopoEntityRef) {
        if let Some(old_target) = self.name_to_ref.remove(&name) {
            self.ref_to_name.remove(&old_target);
        }
        if let Some(old_name) = self.ref_to_name.remove(&target) {
            self.name_to_ref.remove(&old_name);
        }
        self.name_to_ref.insert(name.clone(), target);
        self.ref_to_name.insert(target, name);
    }
}

fn flat_face_count(brep: &BRep) -> usize {
    brep
        .solids
        .iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum()
}

fn is_valid_ref_for_brep(brep: &BRep, target: TopoEntityRef) -> bool {
    match target {
        TopoEntityRef::Vertex(i) => i < brep.vertices.len(),
        TopoEntityRef::Edge(i) => i < brep.edges.len(),
        TopoEntityRef::Face(i) => i < flat_face_count(brep),
        TopoEntityRef::Solid(i) => i < brep.solids.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BRep, PrimitiveSolid};

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    #[test]
    fn default_labels_cover_basic_topology() {
        let brep = unit_box();
        let hooks = PersistentNamingHooks::with_default_labels_for_brep(&brep);
        assert_eq!(hooks.resolve("v0"), Some(TopoEntityRef::Vertex(0)));
        assert_eq!(hooks.resolve("e0"), Some(TopoEntityRef::Edge(0)));
        assert_eq!(hooks.resolve("f0"), Some(TopoEntityRef::Face(0)));
        assert_eq!(hooks.resolve("s0"), Some(TopoEntityRef::Solid(0)));
    }

    #[test]
    fn bind_and_rename_roundtrip() {
        let brep = unit_box();
        let mut hooks = PersistentNamingHooks::new();
        hooks
            .bind_for_brep(&brep, "mount_hole", TopoEntityRef::Edge(1))
            .expect("bind should succeed");
        assert_eq!(hooks.resolve("mount_hole"), Some(TopoEntityRef::Edge(1)));

        hooks
            .rename("mount_hole", "outer_profile")
            .expect("rename should succeed");
        assert_eq!(hooks.resolve("mount_hole"), None);
        assert_eq!(hooks.resolve("outer_profile"), Some(TopoEntityRef::Edge(1)));
    }

    #[test]
    fn validate_and_retain_invalid_bindings() {
        let brep = unit_box();
        let mut hooks = PersistentNamingHooks::new();
        hooks.bind("bad_edge", TopoEntityRef::Edge(9999));
        hooks.bind("good_vertex", TopoEntityRef::Vertex(0));

        let issues = hooks.validate_against_brep(&brep);
        assert_eq!(issues.len(), 1);
        assert!(issues[0].contains("bad_edge"));

        hooks.retain_valid_for_brep(&brep);
        assert_eq!(hooks.resolve("bad_edge"), None);
        assert_eq!(hooks.resolve("good_vertex"), Some(TopoEntityRef::Vertex(0)));
    }
}
