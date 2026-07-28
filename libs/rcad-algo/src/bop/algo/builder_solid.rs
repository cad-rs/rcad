// OCCT BOPAlgo_BuilderSolid — solid building from shells.
//
// OCCT BOPAlgo_BuilderSolid.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use rcad_kernel::topo_shape::Shape;

/// OCCT BOPAlgo_BuilderSolid — builds solids from a set of faces.
pub struct BuilderSolid<'a> {
    ds: &'a DS,
    my_report: Report,
    /// Input faces (OCCT: myShapes)
    pub my_shapes: Vec<Shape>,
    /// Resulting solids (OCCT: mySolids)
    pub my_solids: Vec<Shape>,
}

impl<'a> BuilderSolid<'a> {
    pub fn new(ds: &'a DS) -> Self {
        BuilderSolid {
            ds,
            my_report: Report::new(),
            my_shapes: Vec::new(),
            my_solids: Vec::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::Perform (BOPAlgo_BuilderSolid.cxx L76-125).
    pub fn perform(&mut self) {
        if self.my_shapes.is_empty() {
            return;
        }
        // OCCT L106: PerformShapesToAvoid
        self.perform_shapes_to_avoid();
        if self.has_errors() { return; }
        // OCCT L112: PerformLoops
        self.perform_loops();
        if self.has_errors() { return; }
        // OCCT L118: PerformAreas
        self.perform_areas();
        if self.has_errors() { return; }
        // OCCT L124: PerformInternalShapes
        self.perform_internal_shapes();
    }

    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }

    /// OCCT BOPAlgo_BuilderSolid::PerformShapesToAvoid (BOPAlgo_BuilderSolid.cxx L129+).
    fn perform_shapes_to_avoid(&mut self) {}

    /// OCCT BOPAlgo_BuilderSolid::PerformLoops.
    fn perform_loops(&mut self) {}

    /// OCCT BOPAlgo_BuilderSolid::PerformAreas.
    fn perform_areas(&mut self) {}

    /// OCCT BOPAlgo_BuilderSolid::PerformInternalShapes.
    fn perform_internal_shapes(&mut self) {}
}
