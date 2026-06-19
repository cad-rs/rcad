//! OCCT-aligned BuilderSolid: builds closed solids from a set of faces.
//!
//! OCCT ref: BOPAlgo_BuilderSolid (BOPAlgo_BuilderSolid.cxx / .hxx)
//!
//! OCCT BOPAlgo_BuilderSolid::Perform steps:
//! 1. PerformShapesToAvoid — identify internal/section faces to exclude
//! 2. PerformLoops — build oriented loops from section edges
//! 3. PerformAreas — classify 3D areas as inside/outside
//! 4. BuildSolids — construct TopoDS_Solid from classified areas
//!
//! ⏳ Partial alignment: currently uses ShellSplitter to find connected
//!   components (step 0).  The full PerformLoops/PerformAreas/BuildSolids
//!   pipeline is a future target.
//!
//! NOTE: In OCCT, BuilderSolid can produce multiple solids from disconnected
//! face sets.  The current rcad boolean builder creates a single solid;
//! this struct bridges the gap by splitting disconnected components.

use super::ds::DS;
use super::shell_splitter::ShellSplitter;

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct BuilderSolid {
    /// Shell splitters, one per connected component of faces.
    /// OCCT: BOPAlgo_BuilderSolid creates one ShellSplitter per block.
    myShells: Vec<ShellSplitter>,
    /// Tolerance used for geometric classification.
    /// OCCT: inherited from BOPAlgo_Options (BOPAlgo_BuilderSolid inherits
    /// BOPAlgo_Builder via BOPAlgo_BuilderShape).
    myTolerance: f64,
    /// Resulting solids: each entry is a Vec of face indices forming a closed solid.
    /// OCCT: BOPAlgo_BuilderSolid::mySolids (list of TopoDS_Solid).
    mySolids: Vec<Vec<usize>>,
}

impl BuilderSolid {
    /// Create an empty BuilderSolid.
    ///
    /// ✅ OCCT-aligned: BOPAlgo_BuilderSolid default constructor.
    pub fn new() -> Self {
        Self {
            myShells: Vec::new(),
            myTolerance: 1e-7,
            mySolids: Vec::new(),
        }
    }

    /// Set the tolerance for geometric classification.
    ///
    /// ✅ OCCT-aligned: BOPAlgo_Options::SetTolerance.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.myTolerance = tol;
    }

    /// Build closed solids from a set of faces.
    ///
    /// ✅ OCCT-aligned: BOPAlgo_BuilderSolid::Perform (top-level dispatch,
    ///   BOPAlgo_BuilderSolid.cxx lines 1-60).
    ///
    /// Current implementation:
    ///   1. Uses ShellSplitter to find connected face components.
    ///   2. Each connected component becomes a solid candidate.
    ///
    /// Future (full OCCT alignment):
    ///   - PerformShapesToAvoid: filter internal faces (faces whose normals
    ///     point inward or faces from intersection history).
    ///   - PerformLoops: build oriented wire loops from section edges on
    ///     each shell (BOPAlgo_BuilderSolid_Loops.cxx).
    ///   - PerformAreas: classify bounded 3D regions as in/out using
    ///     BOPAlgo_Tools::ClassifyFaces (BOPAlgo_BuilderSolid_Areas.cxx).
    ///   - BuildSolids: construct TopoDS_Solid per classified region
    ///     (BOPAlgo_BuilderSolid_Solid.cxx).
    pub fn perform(&mut self, ds: &DS, all_faces: &[usize]) {
        self.mySolids.clear();
        self.myShells.clear();

        if all_faces.is_empty() {
            return;
        }

        // Step 1: Partition faces into connected components (shells).
        // OCCT ref: BOPAlgo_ShellSplitter is used as a building block
        // (BOPAlgo_BuilderSolid creates a ShellSplitter internally).
        let mut splitter = ShellSplitter::new();
        for &fi in all_faces {
            splitter.add_start_face(fi);
        }
        splitter.perform(ds);

        // Step 2 (stub): Each connected component becomes a solid candidate.
        // OCCT ref: PerformLoops + PerformAreas would classify each shell's
        // interior/exterior before building solids.
        for shell_faces in splitter.shells() {
            if !shell_faces.is_empty() {
                self.mySolids.push(shell_faces.clone());
            }
        }

        self.myShells.push(splitter);
    }

    /// Return the resulting solids (each entry is a vector of face indices).
    ///
    /// ✅ OCCT-aligned: BOPAlgo_BuilderSolid::Solids() accessor.
    pub fn solids(&self) -> &[Vec<usize>] {
        &self.mySolids
    }

    /// Return the number of solids built.
    pub fn nb_solids(&self) -> usize {
        self.mySolids.len()
    }

    /// Return true when multiple solids were produced.
    pub fn has_multiple_solids(&self) -> bool {
        self.mySolids.len() > 1
    }

    /// Access the shell splitters used internally.
    pub fn shell_splitters(&self) -> &[ShellSplitter] {
        &self.myShells
    }
}

impl Default for BuilderSolid {
    fn default() -> Self {
        Self::new()
    }
}
