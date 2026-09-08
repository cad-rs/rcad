//! OCCT BRepSweep_Tool (TKPrim/BRepSweep) — the indexation and type analysis
//! services required by the TopoDS generating Shape of BRepSweep.
//!
//! Sources:
//! - BRepSweep_Tool.hxx L32-64
//! - BRepSweep_Tool.cxx L23-72
//!
//! Architecture difference: OCCT keeps an NCollection_IndexedMap of the
//! shapes (TopExp::MapShapes over TShape+Location identity); the rcad map is
//! an IndexMap keyed by the TopTools_ShapeMapHasher pair (TShape pointer +
//! Location, orientation ignored).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use indexmap::IndexMap;

use crate::brep_algo::tool::{shape_key, sub_shapes};

pub(crate) type ShapeKey = (u64, u32);

/// OCCT BRepSweep_Tool (BRepSweep_Tool.hxx L32-62).
pub struct BRepSweepTool {
    my_map: IndexMap<ShapeKey, Shape>,
}

impl BRepSweepTool {
    /// OCCT BRepSweep_Tool::BRepSweep_Tool(aShape) (cxx L23-26) —
    /// TopExp::MapShapes(aShape, myMap): the shape and all its sub-shapes,
    /// outer first, inserted in the exploration order.
    pub fn new(a_shape: &Shape) -> Self {
        let mut my_map = IndexMap::new();
        map_shapes(a_shape, &mut my_map);
        BRepSweepTool { my_map }
    }

    /// OCCT BRepSweep_Tool::NbShapes() (cxx L30-33) — the number of subshapes
    /// in the shape (myMap.Extent()).
    pub fn nb_shapes(&self) -> i32 {
        self.my_map.len() as i32
    }

    /// OCCT BRepSweep_Tool::Index(aShape) (cxx L37-44) — the index of
    /// <aShape> (0 when not contained).
    pub fn index(&self, a_shape: &Shape) -> i32 {
        match self.my_map.get_index_of(&shape_key(a_shape)) {
            Some(i) => (i + 1) as i32,
            None => 0,
        }
    }

    /// OCCT BRepSweep_Tool::Shape(anIndex) (cxx L48-51) — the Shape at Index
    /// anIndex (myMap.FindKey(anIndex)).
    pub fn shape(&self, an_index: i32) -> Shape {
        self.my_map
            .get_index((an_index - 1) as usize)
            .expect("myMap.FindKey")
            .1
            .clone()
    }

    /// OCCT BRepSweep_Tool::Type(aShape) (cxx L55-58) — the type of <aShape>.
    pub fn type_of(&self, a_shape: &Shape) -> ShapeType {
        a_shape.shape_type()
    }

    /// OCCT BRepSweep_Tool::Orientation(aShape) (cxx L62-65) — the
    /// Orientation of <aShape>.
    pub fn orientation(&self, a_shape: &Shape) -> Orientation {
        a_shape.orientation
    }

    /// OCCT BRepSweep_Tool::SetOrientation(aShape, Or) (cxx L69-72) — sets
    /// the Orientation of <aShape> with Or (the handle's orientation field,
    /// no TShape change).
    pub fn set_orientation(&self, a_shape: &mut Shape, or: Orientation) {
        a_shape.orientation = or;
    }
}

/// OCCT TopExp::MapShapes(S, M) (TopExp.cxx MapShapes): M.Add(S) then the
/// recursive walk over the direct sub-shapes (TopoDS_Iterator order).
fn map_shapes(s: &Shape, m: &mut IndexMap<ShapeKey, Shape>) {
    if m.insert(shape_key(s), s.clone()).is_some() {
        return;
    }
    for c in sub_shapes(s) {
        map_shapes(&c, m);
    }
}

