//! OCCT Sweep_NumShapeTool (TKPrim/Sweep) — the indexation and type analysis
//! services required by the NumShape Directing Shapes of Swept Primitives.
//!
//! Sources:
//! - Sweep_NumShapeTool.hxx L31-72
//! - Sweep_NumShapeTool.cxx L23-154

use rcad_kernel::topods::{Orientation, ShapeType};

use super::sweep_num_shape::SweepNumShape;

/// OCCT Sweep_NumShapeTool (Sweep_NumShapeTool.hxx L31-70).
#[derive(Debug, Clone, Copy)]
pub struct SweepNumShapeTool {
    my_num_shape: SweepNumShape,
}

impl SweepNumShapeTool {
    /// OCCT Sweep_NumShapeTool::Sweep_NumShapeTool(aShape) (cxx L23-26) —
    /// creates a new NumShapeTool with <aShape>.
    pub fn new(a_shape: &SweepNumShape) -> Self {
        SweepNumShapeTool {
            my_num_shape: *a_shape,
        }
    }

    /// OCCT Sweep_NumShapeTool::NbShapes() (cxx L30-47) — the number of
    /// subshapes in the shape.
    pub fn nb_shapes(&self) -> i32 {
        if self.my_num_shape.type_() == ShapeType::Edge {
            if self.my_num_shape.closed() {
                self.my_num_shape.index()
            } else {
                self.my_num_shape.index() + 1
            }
        } else {
            1
        }
    }

    /// OCCT Sweep_NumShapeTool::Index(aShape) (cxx L51-68) — the index of
    /// <aShape>.
    pub fn index(&self, a_shape: &SweepNumShape) -> i32 {
        if a_shape.type_() == ShapeType::Edge {
            1
        } else if a_shape.closed() {
            2
        } else {
            a_shape.index() + 1
        }
    }

    /// OCCT Sweep_NumShapeTool::Shape(anIndex) (cxx L72-82) — the Shape at
    /// index anIndex.
    pub fn shape(&self, an_index: i32) -> SweepNumShape {
        if an_index == 1 {
            self.my_num_shape
        } else {
            SweepNumShape::with_all(
                an_index - 1,
                ShapeType::Vertex,
                self.my_num_shape.closed(),
                false,
                false,
            )
        }
    }

    /// OCCT Sweep_NumShapeTool::Type(aShape) (cxx L86-89) — the type of
    /// <aShape>.
    pub fn type_of(&self, a_shape: &SweepNumShape) -> ShapeType {
        a_shape.type_()
    }

    /// OCCT Sweep_NumShapeTool::Orientation(aShape) (cxx L93-96) — the
    /// orientation of <aShape>.
    pub fn orientation(&self, a_shape: &SweepNumShape) -> Orientation {
        a_shape.orientation()
    }

    /// OCCT Sweep_NumShapeTool::HasFirstVertex() (cxx L100-107) — true if
    /// there is a First Vertex in the Shape.
    pub fn has_first_vertex(&self) -> bool {
        if self.my_num_shape.type_() == ShapeType::Edge {
            return !self.my_num_shape.beg_infinite();
        }
        true
    }

    /// OCCT Sweep_NumShapeTool::HasLastVertex() (cxx L111-118) — true if
    /// there is a Last Vertex in the Shape.
    pub fn has_last_vertex(&self) -> bool {
        if self.my_num_shape.type_() == ShapeType::Edge {
            return !self.my_num_shape.end_infinite();
        }
        true
    }

    /// OCCT Sweep_NumShapeTool::FirstVertex() (cxx L122-136) — the first
    /// vertex.
    pub fn first_vertex(&self) -> SweepNumShape {
        if self.my_num_shape.type_() == ShapeType::Edge {
            if self.has_first_vertex() {
                return SweepNumShape::with_all(
                    1,
                    ShapeType::Vertex,
                    self.my_num_shape.closed(),
                    false,
                    false,
                );
            } else {
                // OCCT: throw Standard_ConstructionError("infinite Shape").
                panic!("Standard_ConstructionError: infinite Shape");
            }
        }
        self.my_num_shape
    }

    /// OCCT Sweep_NumShapeTool::LastVertex() (cxx L140-154) — the last
    /// vertex.
    pub fn last_vertex(&self) -> SweepNumShape {
        if self.my_num_shape.type_() == ShapeType::Edge {
            if self.has_last_vertex() {
                return SweepNumShape::with_all(
                    self.nb_shapes() - 1,
                    ShapeType::Vertex,
                    self.my_num_shape.closed(),
                    false,
                    false,
                );
            } else {
                // OCCT: throw Standard_ConstructionError("infinite Shape").
                panic!("Standard_ConstructionError: infinite Shape");
            }
        }
        self.my_num_shape
    }
}
