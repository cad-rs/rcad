// BRepTools_History semantic equivalent (BRepTools_History.cxx / .hxx) — the
// history of modifications, generations and removals of the boolean operation.

use rcad_kernel::topods::{ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use std::collections::HashMap;
use std::collections::HashSet;

/// OCCT BRepTools_History::IsSupportedType (BRepTools_History.hxx L145-153):
/// only VERTEX, EDGE, FACE and SOLID are supported.
pub fn is_supported_type(s: &Shape) -> bool {
    match s.shape_type() {
        ShapeType::Vertex | ShapeType::Edge | ShapeType::Face | ShapeType::Solid => true,
        _ => false,
    }
}

/// OCCT BRepTools_History (BRepTools_History.cxx) — myShapeToModified /
/// myShapeToGenerated DataMaps and the myRemoved map.
pub struct BRepToolsHistory {
    // myShapeToModified: DataMap<Shape, List<Shape>>
    my_shape_to_modified: HashMap<(u64, u32), Vec<Shape>>,
    // myShapeToGenerated: DataMap<Shape, List<Shape>>
    my_shape_to_generated: HashMap<(u64, u32), Vec<Shape>>,
    // myRemoved: NCollection_Map<Shape> (TShape + Location)
    my_removed: HashSet<(u64, u32)>,
}

impl BRepToolsHistory {
    pub fn new() -> Self {
        Self {
            my_shape_to_modified: HashMap::new(),
            my_shape_to_generated: HashMap::new(),
            my_removed: HashSet::new(),
        }
    }

    fn key(s: &Shape) -> (u64, u32) {
        (s.ptr_id(), s.location)
    }

    /// OCCT BRepTools_History::AddGenerated (BRepTools_History.cxx L48-67).
    pub fn add_generated(&mut self, the_initial: &Shape, the_generated: &Shape) {
        if !is_supported_type(the_initial) || !is_supported_type(the_generated) {
            return;
        }
        let list = self
            .my_shape_to_generated
            .entry(Self::key(the_initial))
            .or_default();
        if !list.iter().any(|g| Self::key(g) == Self::key(the_generated)) {
            list.push(the_generated.clone());
        }
    }

    /// OCCT BRepTools_History::AddModified (BRepTools_History.cxx L69-88).
    pub fn add_modified(&mut self, the_initial: &Shape, the_modified: &Shape) {
        if !is_supported_type(the_initial) || !is_supported_type(the_modified) {
            return;
        }
        let list = self
            .my_shape_to_modified
            .entry(Self::key(the_initial))
            .or_default();
        if !list.iter().any(|m| Self::key(m) == Self::key(the_modified)) {
            list.push(the_modified.clone());
        }
    }

    /// OCCT BRepTools_History::Remove (BRepTools_History.cxx L91-108) — unbind
    /// the modifications and add the shape to myRemoved.
    pub fn remove(&mut self, the_removed: &Shape) {
        if !is_supported_type(the_removed) {
            return;
        }
        self.my_shape_to_modified.remove(&Self::key(the_removed));
        self.my_removed.insert(Self::key(the_removed));
    }

    /// OCCT BRepTools_History::Modified.
    pub fn modified(&self, the_initial: &Shape) -> Vec<Shape> {
        self.my_shape_to_modified
            .get(&Self::key(the_initial))
            .cloned()
            .unwrap_or_default()
    }

    /// OCCT BRepTools_History::Generated.
    pub fn generated(&self, the_initial: &Shape) -> Vec<Shape> {
        self.my_shape_to_generated
            .get(&Self::key(the_initial))
            .cloned()
            .unwrap_or_default()
    }

    /// OCCT BRepTools_History::IsRemoved.
    pub fn is_removed(&self, the_initial: &Shape) -> bool {
        self.my_removed.contains(&Self::key(the_initial))
    }

    /// Sanity: TShape referenced for dead-code elimination.
    #[allow(dead_code)]
    fn _shape_type(s: &Shape) -> Option<&'static str> {
        match &*s.data {
            TShape::Vertex(_) => Some("vertex"),
            TShape::Edge(_) => Some("edge"),
            TShape::Face(_) => Some("face"),
            TShape::Solid(_) => Some("solid"),
            _ => None,
        }
    }
}
