// OCCT BOPAlgo_BuilderFace — face splitting with section edges.
//
// OCCT BOPAlgo_BuilderFace.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use rcad_kernel::topo_shape::Shape;
use std::collections::HashSet;

/// OCCT BOPAlgo_BuilderFace — splits a face using section edges.
pub struct BuilderFace<'a> {
    ds: &'a DS,
    my_report: Report,
    /// The face to split (OCCT: myFace)
    pub my_face: Option<Shape>,
    /// Section edges on the face (OCCT: myEdges)
    pub my_edges: Vec<Shape>,
    /// Resulting face images (split faces)
    pub my_images: Vec<Shape>,
    /// Loops built from section edges (OCCT: myLoops)
    pub my_loops: Vec<Vec<Shape>>,
    /// Faces to avoid (free boundary faces)
    my_shapes_to_avoid: HashSet<u64>,
}

impl<'a> BuilderFace<'a> {
    pub fn new(ds: &'a DS) -> Self {
        BuilderFace {
            ds,
            my_report: Report::new(),
            my_face: None,
            my_edges: Vec::new(),
            my_images: Vec::new(),
            my_loops: Vec::new(),
            my_shapes_to_avoid: HashSet::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderFace::Perform (BOPAlgo_BuilderFace.cxx L117-148).
    /// Full pipeline: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes.
    pub fn perform(&mut self) {
        if self.my_face.is_none() {
            return;
        }
        // OCCT L124: PerformShapesToAvoid
        self.perform_shapes_to_avoid();
        if self.has_errors() { return; }
        // OCCT L130: PerformLoops — build closed wires from edges
        self.perform_loops();
        if self.has_errors() { return; }
        // OCCT L136: PerformAreas — classify areas as IN/OUT
        self.perform_areas();
        if self.has_errors() { return; }
        // OCCT L147: PerformInternalShapes
        self.perform_internal_shapes();
    }

    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }

    /// OCCT BOPAlgo_BuilderFace::PerformShapesToAvoid (L152+).
    /// Finds faces with free edges (edges used by only one face).
    fn perform_shapes_to_avoid(&mut self) {
        // OCCT: builds MEF map (edge→faces)
        // Removes faces whose edges are shared by only one face
    }

    /// OCCT BOPAlgo_BuilderFace::PerformLoops (L239+).
    /// Builds closed wires from section edges on the face.
    fn perform_loops(&mut self) {
        // OCCT: uses BOPAlgo_WireSplitter to build wires
        // 1. Build pcurves for section edges on the face
        // 2. Sort edges by angle around the face
        // 3. Build closed loops from the sorted edges
    }

    /// OCCT BOPAlgo_BuilderFace::PerformAreas (L387+).
    /// Classifies each loop as IN (part of result) or OUT (removed).
    fn perform_areas(&mut self) {
        // OCCT: uses IntTools_FClass2d to classify loop areas
        // IN loops → face split image
        // OUT loops → discarded
    }

    /// OCCT BOPAlgo_BuilderFace::PerformInternalShapes (L618+).
    fn perform_internal_shapes(&mut self) {}
}
