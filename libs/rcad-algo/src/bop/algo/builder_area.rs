//! OCCT BOPAlgo_BuilderArea — root class for building faces/solids from edges/faces.
//!
//! OCCT BOPAlgo_BuilderArea.cxx / .hxx:
//! - myContext, myShapes, myLoops, myLoopsInternal, myAreas, myShapesToAvoid,
//!   myAvoidInternalShapes.
//! - Accessors SetContext/Shapes/SetShapes/Loops/Areas/SetAvoidInternalShapes/
//!   IsAvoidInternalShapes.
//! - Pure virtual PerformShapesToAvoid/PerformLoops/PerformAreas/PerformInternalShapes.
//!
//! Translation notes:
//! - NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! - NCollection_IndexedMap<TopoDS_Shape> -> HashSet<u64> keyed by ptr_id
//!   (TopTools_ShapeMapHasher = TShape + Location, orientation-insensitive).
//! - Rust has no inheritance: BuilderFace/BuilderSolid hold the same members
//!   directly (my_loops/my_loops_internal as Vec<Vec<Shape>> edge sequences).

use rcad_kernel::topo_shape::Shape;
use crate::bop::algo::Report;
use crate::bop::int_tools::context::IntToolsContext;
use std::collections::HashSet;

/// OCCT BOPAlgo_BuilderArea (BOPAlgo_BuilderArea.hxx).
pub struct BuilderArea {
    /// OCCT: myContext.
    pub my_context: Option<IntToolsContext>,
    /// OCCT: myShapes — the input shapes.
    pub my_shapes: Vec<Shape>,
    /// OCCT: myLoops — the found loops (wires).
    pub my_loops: Vec<Shape>,
    /// OCCT: myLoopsInternal — the internal loops (wires).
    pub my_loops_internal: Vec<Shape>,
    /// OCCT: myAreas — the found areas (result faces/solids).
    pub my_areas: Vec<Shape>,
    /// OCCT: myShapesToAvoid — shapes to be excluded from the loops.
    pub my_shapes_to_avoid: HashSet<u64>,
    /// OCCT: myAvoidInternalShapes — prevents addition of internal parts.
    pub my_avoid_internal_shapes: bool,
    /// BOPAlgo_Algo (inherited).
    pub my_report: Report,
    pub my_run_parallel: bool,
    pub my_fuzzy_value: f64,
}

impl BuilderArea {
    pub fn new() -> Self {
        Self {
            my_context: None,
            my_shapes: Vec::new(),
            my_loops: Vec::new(),
            my_loops_internal: Vec::new(),
            my_areas: Vec::new(),
            my_shapes_to_avoid: HashSet::new(),
            my_avoid_internal_shapes: false,
            my_report: Report::new(),
            my_run_parallel: false,
            // OCCT BOPAlgo_Algo base: myFuzzyValue(Precision::Confusion()).
            my_fuzzy_value: rcad_kernel::precision::CONFUSION,
        }
    }

    /// OCCT: SetContext(theContext).
    pub fn set_context(&mut self, ctx: IntToolsContext) {
        self.my_context = Some(ctx);
    }

    /// OCCT: SetShapes(theLS).
    pub fn set_shapes(&mut self, shapes: Vec<Shape>) {
        self.my_shapes = shapes;
    }

    /// OCCT: Shapes().
    pub fn shapes(&self) -> &[Shape] {
        &self.my_shapes
    }

    /// OCCT: Loops().
    pub fn loops(&self) -> &[Shape] {
        &self.my_loops
    }

    /// OCCT: Areas().
    pub fn areas(&self) -> &[Shape] {
        &self.my_areas
    }

    /// OCCT: SetAvoidInternalShapes(theAvoidInternal).
    pub fn set_avoid_internal_shapes(&mut self, the_avoid: bool) {
        self.my_avoid_internal_shapes = the_avoid;
    }

    /// OCCT: IsAvoidInternalShapes().
    pub fn is_avoid_internal_shapes(&self) -> bool {
        self.my_avoid_internal_shapes
    }
}
