use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use crate::topo::topods;
use crate::topo::topo_query::vertex_indices;

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
    /// - semantic vertices: `v0`, `v1`, ...
    /// - edges: `e0`, `e1`, ...
    /// - faces: `f0`, `f1`, ... (tshape index)
    /// - solids: `s0`, `s1`, ...
    pub fn with_default_labels_for_brep(brep: &topods::BRep) -> Self {
        let mut out = Self::new();
        for vi in vertex_indices(brep) {
            out.bind_unchecked(format!("v{vi}"), TopoEntityRef::Vertex(vi));
        }
        for (i, ts) in brep.tshapes.iter().enumerate() {
            match &**ts {
                topods::TShape::Edge(_) => {
                    out.bind_unchecked(format!("e{i}"), TopoEntityRef::Edge(i));
                }
                topods::TShape::Face(_) => {
                    out.bind_unchecked(format!("f{i}"), TopoEntityRef::Face(i));
                }
                topods::TShape::Solid(_) => {
                    out.bind_unchecked(format!("s{i}"), TopoEntityRef::Solid(i));
                }
                _ => {}
            }
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
        brep: &topods::BRep,
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
    pub fn validate_against_brep(&self, brep: &topods::BRep) -> Vec<String> {
        let mut issues = Vec::new();
        for (name, target) in &self.name_to_ref {
            if !is_valid_ref_for_brep(brep, *target) {
                issues.push(format!(
                    "name '{name}' points to out-of-range entity {target:?}"
                ));
            }
        }
        issues
    }

    /// Remove bindings that no longer point to valid topology entities.
    pub fn retain_valid_for_brep(&mut self, brep: &topods::BRep) {
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

    /// Propagate names from a pre-operation naming table through an index-level
    /// mapping to produce a post-operation naming table.
    ///
    /// `face_map[old_face_idx]` 鈫?`Some(new_face_idx)` if the face survived the
    /// operation and was remapped.  `None` means the face was consumed/deleted.
    /// New entities (generated by the operation) can be named via the
    /// `new_face_names` slice: `(new_face_idx, name)` pairs.
    ///
    /// The same logic applies independently to vertices and edges via their
    /// respective maps.  This method updates `self` in place and returns the set
    /// of names that were dropped (because the entity was removed).
    ///
    /// Analogous to OCCT `TNaming_NamedShape` propagation after a BRep rebuild.
    pub fn propagate_through_remap(
        &mut self,
        face_map: &[Option<usize>],
        edge_map: &[Option<usize>],
        vertex_map: &[Option<usize>],
        new_face_names: &[(usize, String)],
        new_edge_names: &[(usize, String)],
        new_vertex_names: &[(usize, String)],
    ) -> Vec<String> {
        let mut dropped = Vec::new();

        // Collect existing bindings so we can update them.
        let snapshot: Vec<(String, TopoEntityRef)> = self
            .name_to_ref
            .iter()
            .map(|(n, r)| (n.clone(), *r))
            .collect();

        // Clear tables; we rebuild them below.
        self.name_to_ref.clear();
        self.ref_to_name.clear();

        for (name, old_ref) in snapshot {
            let new_ref = match old_ref {
                TopoEntityRef::Face(i) => {
                    if face_map.is_empty() {
                        Some(old_ref)
                    } else {
                        face_map.get(i).and_then(|r| *r).map(TopoEntityRef::Face)
                    }
                }
                TopoEntityRef::Edge(i) => {
                    if edge_map.is_empty() {
                        Some(old_ref)
                    } else {
                        edge_map.get(i).and_then(|r| *r).map(TopoEntityRef::Edge)
                    }
                }
                TopoEntityRef::Vertex(i) => {
                    if vertex_map.is_empty() {
                        Some(old_ref)
                    } else {
                        vertex_map
                            .get(i)
                            .and_then(|r| *r)
                            .map(TopoEntityRef::Vertex)
                    }
                }
                TopoEntityRef::Solid(_) => Some(old_ref), // solids not remapped
            };
            match new_ref {
                Some(r) => self.bind_unchecked(name, r),
                None => dropped.push(name),
            }
        }

        // Register names for new entities.
        for (idx, name) in new_face_names {
            self.bind_unchecked(name.clone(), TopoEntityRef::Face(*idx));
        }
        for (idx, name) in new_edge_names {
            self.bind_unchecked(name.clone(), TopoEntityRef::Edge(*idx));
        }
        for (idx, name) in new_vertex_names {
            self.bind_unchecked(name.clone(), TopoEntityRef::Vertex(*idx));
        }

        dropped
    }

    /// Convenience: propagate names from an operation that only remaps faces
    /// (e.g. fillet, chamfer).
    ///
    /// `face_map[old_idx]` 鈫?`Some(new_idx)` or `None` (removed).
    /// `new_face_names` names any newly generated faces.
    pub fn propagate_face_remap(
        &mut self,
        face_map: &[Option<usize>],
        new_face_names: &[(usize, String)],
    ) -> Vec<String> {
        self.propagate_through_remap(face_map, &[], &[], new_face_names, &[], &[])
    }

    /// Build a simple identity remap for `n` entities (nothing moved or removed).
    pub fn identity_map(n: usize) -> Vec<Option<usize>> {
        (0..n).map(Some).collect()
    }

    /// Iterate over all (name, entity_ref) bindings.
    pub fn iter(&self) -> impl Iterator<Item = (&str, TopoEntityRef)> {
        self.name_to_ref.iter().map(|(n, r)| (n.as_str(), *r))
    }
}

fn flat_face_count(brep: &topods::BRep) -> usize {
    brep.tshapes
        .iter()
        .filter(|ts| matches!((&**ts).as_ref(), &topods::TShape::Face(_)))
        .count()
}

fn is_valid_ref_for_brep(brep: &topods::BRep, target: TopoEntityRef) -> bool {
    match target {
        TopoEntityRef::Vertex(i) => brep
            .tshapes
            .get(i)
            .is_some_and(|ts| matches!(&**ts, &topods::TShape::Vertex(_))),
        TopoEntityRef::Edge(i) => brep
            .tshapes
            .get(i)
            .is_some_and(|ts| matches!(&**ts, topods::TShape::Edge(_))),
        TopoEntityRef::Face(i) => brep
            .tshapes
            .get(i)
            .is_some_and(|ts| matches!(&**ts, topods::TShape::Face(_))),
        TopoEntityRef::Solid(i) => brep
            .tshapes
            .get(i)
            .is_some_and(|ts| matches!(&**ts, topods::TShape::Solid(_))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn brep_with_edge() -> (topods::BRep, usize) {
        let mut brep = topods::BRep::new();
        // Vertices must come first in tshapes (add_tvertex pushes).
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = brep.add_tedge(
            None,
            v0,
            topods::Shape::synthetic(v1.index, topods::Orientation::Reversed),
            [0.0, 1.0],
        );
        // Edge is at tshapes[2] (after 2 vertices). TopoEntityRef::Edge(i) uses tshapes index.
        (brep, e.index)
    }

    #[test]
    fn bind_and_rename_roundtrip() {
        let mut hooks = PersistentNamingHooks::new();
        let (brep, edge_ts_idx) = brep_with_edge();
        hooks
            .bind_for_brep(&brep, "test", TopoEntityRef::Edge(edge_ts_idx))
            .expect("bind should succeed");
        assert_eq!(
            hooks.resolve("test"),
            Some(TopoEntityRef::Edge(edge_ts_idx))
        );

        hooks
            .rename("test", "renamed")
            .expect("rename should succeed");
        assert_eq!(hooks.resolve("test"), None);
        assert_eq!(
            hooks.resolve("renamed"),
            Some(TopoEntityRef::Edge(edge_ts_idx))
        );
    }

    #[test]
    fn invalid_ref_rejected() {
        let brep = topods::BRep::new();
        let mut hooks = PersistentNamingHooks::new();
        // Edge(0) is invalid because BRep has no edges
        assert!(
            hooks
                .bind_for_brep(&brep, "bad_edge", TopoEntityRef::Edge(0))
                .is_err()
        );
    }

    #[test]
    fn default_labels_empty_brep() {
        let brep = topods::BRep::new();
        let hooks = PersistentNamingHooks::with_default_labels_for_brep(&brep);
        assert!(hooks.resolve("v0").is_none());
        assert!(hooks.resolve("e0").is_none());
    }
}
