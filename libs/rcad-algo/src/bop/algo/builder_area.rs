//! OCCT BOPAlgo_BuilderArea — root class for building faces/solids from edges/faces.
//! Stub: inherits Algo, provides SetShapes/SetContext/Shapes/Loops.

use rcad_kernel::topo_shape::Shape;
use crate::bop::algo::Report;
use crate::bop::int_tools::context::IntToolsContext;

/// OCCT BOPAlgo_BuilderArea (BOPAlgo_BuilderArea.hxx).
pub struct BuilderArea {
    pub my_report: Report,
    pub my_run_parallel: bool,
    pub my_fuzzy_value: f64,
    pub my_context: Option<IntToolsContext>,
    pub my_shapes: Vec<Shape>,
    pub my_loops: Vec<Shape>,
}

impl BuilderArea {
    pub fn new() -> Self {
        Self {
            my_report: Report::new(), my_run_parallel: false, my_fuzzy_value: 0.0,
            my_context: None, my_shapes: Vec::new(), my_loops: Vec::new(),
        }
    }
    pub fn set_context(&mut self, ctx: IntToolsContext) { self.my_context = Some(ctx); }
    pub fn set_shapes(&mut self, shapes: Vec<Shape>) { self.my_shapes = shapes; }
    pub fn shapes(&self) -> &[Shape] { &self.my_shapes }
    pub fn loops(&self) -> &[Shape] { &self.my_loops }
}
