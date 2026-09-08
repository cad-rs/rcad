//! OCCT Sweep_NumShapeIterator (TKPrim/Sweep) — iteration services required
//! by the Swept Primitives for a Directing NumShape Line.
//!
//! Sources:
//! - Sweep_NumShapeIterator.hxx L30-62
//! - Sweep_NumShapeIterator.cxx L22-77
//! - Sweep_NumShapeIterator.lxx L19-36 (the inline accessors)

use rcad_kernel::topods::{Orientation, ShapeType};

use super::sweep_num_shape::SweepNumShape;

/// OCCT Sweep_NumShapeIterator (Sweep_NumShapeIterator.hxx L30-58).
#[derive(Debug, Clone, Copy)]
pub struct SweepNumShapeIterator {
    my_num_shape: SweepNumShape,
    my_current_num_shape: SweepNumShape,
    my_current_range: i32,
    my_more: bool,
    my_current_orientation: Orientation,
}

impl Default for SweepNumShapeIterator {
    /// OCCT Sweep_NumShapeIterator::Sweep_NumShapeIterator() (cxx L22-28).
    fn default() -> Self {
        Self::new()
    }
}

impl SweepNumShapeIterator {
    /// OCCT Sweep_NumShapeIterator::Sweep_NumShapeIterator() (cxx L22-28).
    pub fn new() -> Self {
        SweepNumShapeIterator {
            my_num_shape: SweepNumShape::new_indexed(0, ShapeType::Shape),
            my_current_num_shape: SweepNumShape::new_indexed(0, ShapeType::Shape),
            my_current_range: 0,
            my_more: false,
            my_current_orientation: Orientation::Forward,
        }
    }

    /// OCCT Sweep_NumShapeIterator::Init(aShape) (cxx L32-60) — resets the
    /// iterator on the sub-shapes of <aShape>.
    pub fn init(&mut self, a_shape: &SweepNumShape) {
        self.my_num_shape = *a_shape;
        if self.my_num_shape.type_() == ShapeType::Edge {
            let nbvert = self.my_num_shape.index();
            self.my_more = nbvert >= 1;
            if self.my_more {
                self.my_current_range = 1;
                self.my_current_num_shape = SweepNumShape::with_all(
                    1,
                    ShapeType::Vertex,
                    self.my_num_shape.closed(),
                    false,
                    false,
                );
                if nbvert == 1 {
                    if self.my_num_shape.beg_infinite() {
                        self.my_current_orientation = Orientation::Reversed;
                    } else {
                        self.my_current_orientation = Orientation::Forward;
                    }
                } else {
                    self.my_current_orientation = Orientation::Forward;
                }
            }
        }
    }

    /// OCCT Sweep_NumShapeIterator::More() (lxx L19-22) — true if there is a
    /// current sub-shape.
    pub fn more(&self) -> bool {
        self.my_more
    }

    /// OCCT Sweep_NumShapeIterator::Next() (cxx L64-77) — moves to the next
    /// sub-shape.
    pub fn next(&mut self) {
        self.my_current_range += 1;
        self.my_more = self.my_current_range <= self.my_num_shape.index();
        if self.my_more {
            if self.my_num_shape.type_() == ShapeType::Edge {
                self.my_current_num_shape = SweepNumShape::with_all(
                    self.my_current_range,
                    ShapeType::Vertex,
                    self.my_num_shape.closed(),
                    false,
                    false,
                );
                self.my_current_orientation = Orientation::Reversed;
            }
        }
    }

    /// OCCT Sweep_NumShapeIterator::Value() (lxx L26-29) — the current
    /// sub-shape.
    pub fn value(&self) -> &SweepNumShape {
        &self.my_current_num_shape
    }

    /// OCCT Sweep_NumShapeIterator::Orientation() (lxx L33-36) — the
    /// orientation of the current sub-shape.
    pub fn orientation(&self) -> Orientation {
        self.my_current_orientation
    }
}
