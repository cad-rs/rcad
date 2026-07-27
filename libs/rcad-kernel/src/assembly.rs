use glam::DAffine3;
use serde::{Deserialize, Serialize};

use crate::BRep;
use crate::appearance::Color;

/// Assembly shape reference: either a standalone BRep solid or a nested sub-assembly.
///
/// Analogous to an OCCT `XCAFDoc_ShapeTool` shape reference.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(clippy::large_enum_variant)] // BRep leaf is the common case; boxing would churn the API.
pub enum ShapeRef {
    /// Leaf: a single geometric body.
    Brep(BRep),
    /// Non-leaf: nested sub-assembly (`Box` avoids infinite recursive type size).
    Assembly(Box<Assembly>),
}

/// One component instance inside an assembly.
///
/// The same `ShapeRef` may be referenced by multiple `Component`s (instancing);
/// each instance has its own transform and optional color override.
///
/// Analogous to OCCT `NEXT_ASSEMBLY_USAGE_OCCURENCE` + `XCAFDoc_Location`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Component {
    /// Component name (selection, export labels).
    pub name: String,
    /// Geometry content (BRep body or sub-assembly).
    pub shape: ShapeRef,
    /// Transform in the parent assembly frame (translation + rotation + scale).
    pub transform: SerializableAffine3,
    /// Optional color override (overrides the shape’s own color).
    pub color: Option<Color>,
    /// Visibility flag (render filtering).
    pub visible: bool,
}

impl Component {
    /// Create a BRep component at the origin with identity rotation.
    pub fn from_brep(name: impl Into<String>, brep: BRep) -> Self {
        Self {
            name: name.into(),
            shape: ShapeRef::Brep(brep),
            transform: SerializableAffine3::identity(),
            color: None,
            visible: true,
        }
    }

    /// Create a sub-assembly component at the origin with identity rotation.
    pub fn from_assembly(name: impl Into<String>, asm: Assembly) -> Self {
        Self {
            name: name.into(),
            shape: ShapeRef::Assembly(Box::new(asm)),
            transform: SerializableAffine3::identity(),
            color: None,
            visible: true,
        }
    }

    /// Set the transform (builder style).
    pub fn with_transform(mut self, transform: DAffine3) -> Self {
        self.transform = SerializableAffine3(transform);
        self
    }

    /// Set color override (builder style).
    pub fn with_color(mut self, color: Color) -> Self {
        self.color = Some(color);
        self
    }

    /// World-space affine transform for this component.
    pub fn affine(&self) -> DAffine3 {
        self.transform.0
    }
}

/// Hierarchical assembly.
///
/// Holds multiple `Component`s, each with its own transform.
/// Nesting is supported (assemblies inside assemblies).
///
/// Analogous to the shape hierarchy managed by OCCT `XCAFDoc_ShapeTool`.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Assembly {
    /// Assembly name.
    pub name: String,
    /// Child components in order (rendering and traversal).
    pub components: Vec<Component>,
}

impl Assembly {
    /// Create an empty assembly.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            components: Vec::new(),
        }
    }

    /// Append a component.
    pub fn add(&mut self, component: Component) -> &mut Self {
        self.components.push(component);
        self
    }

    /// All BRep leaves in this assembly (flattened) with world transforms.
    ///
    /// Used for rendering, export, collision queries, and any traversal over solid bodies.
    pub fn flatten(&self) -> Vec<FlatComponent> {
        let mut result = Vec::new();
        self.flatten_recursive(DAffine3::IDENTITY, &mut result);
        result
    }

    fn flatten_recursive(&self, parent_transform: DAffine3, out: &mut Vec<FlatComponent>) {
        for component in &self.components {
            if !component.visible {
                continue;
            }
            let world_transform = parent_transform * component.affine();
            match &component.shape {
                ShapeRef::Brep(brep) => {
                    out.push(FlatComponent {
                        name: component.name.clone(),
                        brep: brep.clone(),
                        world_transform,
                        color: component.color,
                    });
                }
                ShapeRef::Assembly(sub_asm) => {
                    sub_asm.flatten_recursive(world_transform, out);
                }
            }
        }
    }

    /// Find a direct or nested component by name (depth-first).
    pub fn find_component(&self, name: &str) -> Option<&Component> {
        for component in &self.components {
            if component.name == name {
                return Some(component);
            }
            if let ShapeRef::Assembly(sub) = &component.shape
                && let Some(found) = sub.find_component(name)
            {
                return Some(found);
            }
        }
        None
    }

    /// Number of immediate child components.
    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// Total BRep leaf count (recursive).
    pub fn leaf_count(&self) -> usize {
        self.components
            .iter()
            .map(|c| match &c.shape {
                ShapeRef::Brep(_) => 1,
                ShapeRef::Assembly(sub) => sub.leaf_count(),
            })
            .sum()
    }
}

/// Flattened leaf: one BRep with a world transform.
#[derive(Debug, Clone)]
pub struct FlatComponent {
    pub name: String,
    pub brep: BRep,
    pub world_transform: DAffine3,
    pub color: Option<Color>,
}

/// Serde-friendly wrapper around `DAffine3` (glam’s type does not implement serde directly).
#[derive(Debug, Clone, Copy)]
pub struct SerializableAffine3(pub DAffine3);

impl SerializableAffine3 {
    pub fn identity() -> Self {
        Self(DAffine3::IDENTITY)
    }
}

impl Serialize for SerializableAffine3 {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        let m = self.0.matrix3;
        let t = self.0.translation;
        // Serialize as 12 f64: [col0 x3, col1 x3, col2 x3, translation x3]
        let mut seq = serializer.serialize_seq(Some(12))?;
        for v in [
            m.x_axis.x, m.x_axis.y, m.x_axis.z, m.y_axis.x, m.y_axis.y, m.y_axis.z, m.z_axis.x,
            m.z_axis.y, m.z_axis.z, t.x, t.y, t.z,
        ] {
            seq.serialize_element(&v)?;
        }
        seq.end()
    }
}

impl<'de> Deserialize<'de> for SerializableAffine3 {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let vals: Vec<f64> = Vec::deserialize(deserializer)?;
        if vals.len() != 12 {
            return Err(serde::de::Error::custom("expected 12 f64 for DAffine3"));
        }
        use glam::{DMat3, DVec3};
        let mat = DMat3::from_cols(
            DVec3::new(vals[0], vals[1], vals[2]),
            DVec3::new(vals[3], vals[4], vals[5]),
            DVec3::new(vals[6], vals[7], vals[8]),
        );
        let trans = DVec3::new(vals[9], vals[10], vals[11]);
        Ok(Self(DAffine3::from_mat3_translation(mat, trans)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BRep;

    #[test]
    fn assembly_flatten() {
        let (box1, _) = BRep::build_unit_cube();
        let (box2, _) = BRep::build_unit_cube();

        let mut asm = Assembly::new("root");
        asm.add(Component::from_brep("box", box1));
        asm.add(Component::from_brep("sphere", box2));

        let flat = asm.flatten();
        assert_eq!(flat.len(), 2);
        assert_eq!(flat[0].name, "box");
        assert_eq!(flat[1].name, "sphere");
    }

    #[test]
    fn assembly_nested() {
        let (box1, _) = BRep::build_unit_cube();

        let mut sub = Assembly::new("sub");
        sub.add(Component::from_brep("part", box1));

        let mut root = Assembly::new("root");
        root.add(Component::from_assembly("sub_asm", sub));

        assert_eq!(root.leaf_count(), 1);
        let flat = root.flatten();
        assert_eq!(flat.len(), 1);
        assert_eq!(flat[0].name, "part");
    }

    #[test]
    fn assembly_find_component() {
        let (box1, _) = BRep::build_unit_cube();
        let mut asm = Assembly::new("root");
        asm.add(Component::from_brep("target", box1));

        assert!(asm.find_component("target").is_some());
        assert!(asm.find_component("missing").is_none());
    }
}
