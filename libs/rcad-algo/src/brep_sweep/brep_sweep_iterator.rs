//! OCCT BRepSweep_Iterator (TKPrim/BRepSweep) — iteration services required
//! by the Generating Line (TopoDS Shape) of a BRepSweep.
//!
//! Sources:
//! - BRepSweep_Iterator.hxx L33-61
//! - BRepSweep_Iterator.cxx L22-36
//! - BRepSweep_Iterator.lxx L19-36 (the inline accessors)
//!
//! Architecture difference: OCCT wraps a TopoDS_Iterator over the
//! TopoDS_Shape children; the rcad walk materializes the children (the
//! TopoDS_Iterator cumOri form, brep_algo::tool::sub_shapes) into a vector
//! with a cursor.

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::Orientation;

/// OCCT BRepSweep_Iterator (BRepSweep_Iterator.hxx L33-57).
pub struct BRepSweepIterator {
    children: Vec<Shape>,
    cursor: usize,
}

impl Default for BRepSweepIterator {
    /// OCCT BRepSweep_Iterator::BRepSweep_Iterator() (cxx L22) — the default
    /// TopoDS_Iterator holds no shape (More() = false).
    fn default() -> Self {
        Self::new()
    }
}

impl BRepSweepIterator {
    /// OCCT BRepSweep_Iterator::BRepSweep_Iterator() (cxx L22) — `= default`.
    pub fn new() -> Self {
        BRepSweepIterator {
            children: Vec::new(),
            cursor: 0,
        }
    }

    /// OCCT BRepSweep_Iterator::Init(aShape) (cxx L26-29) — resets the
    /// iterator on the direct sub-shapes of <aShape>
    /// (myIterator.Initialize(aShape); the TopoDS_Iterator default walks the
    /// direct children with cumulated orientation).
    pub fn init(&mut self, a_shape: &Shape) {
        self.children = super::tool_rehost::topods_sub_shapes(a_shape);
        self.cursor = 0;
    }

    /// OCCT BRepSweep_Iterator::More() (lxx L19-22) — true if there is a
    /// current sub-shape.
    pub fn more(&self) -> bool {
        self.cursor < self.children.len()
    }

    /// OCCT BRepSweep_Iterator::Next() (cxx L33-36) — moves to the next
    /// sub-shape.
    pub fn next(&mut self) {
        self.cursor += 1;
    }

    /// OCCT BRepSweep_Iterator::Value() (lxx L26-29) — the current sub-shape.
    pub fn value(&self) -> Shape {
        self.children[self.cursor].clone()
    }

    /// OCCT BRepSweep_Iterator::Orientation() (lxx L33-36) — the orientation
    /// of the current sub-shape (Value().Orientation()).
    pub fn orientation(&self) -> Orientation {
        self.children[self.cursor].orientation
    }
}
