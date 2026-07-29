//! OCCT BRepAlgoAPI_BuilderOperation equivalent — boolean operation wrapper.
//!
//! Provides a general-purpose boolean operation struct analogous to
//! OCCT BRepAlgoAPI_BuilderOperation (base class for BRepAlgoAPI_Common/Fuse/Cut).
//!
//! OCCT references:
//! - BRepAlgoAPI_Common / BRepAlgoAPI_Fuse / BRepAlgoAPI_Cut (BRepAlgoAPI.cxx)
//! - BRepAlgoAPI_BuilderOperation -> BRepAlgoAPI_BuilderShape -> BRepAlgoAPI_Algo

use crate::bop::algo::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::ds::DS;
use rcad_kernel::topods;
use rcad_kernel::topods::{Orientation, TShape};
use rcad_kernel::topo_shape::Shape;
use std::sync::Arc;

/// A reference to a sub-shape in a rcad_kernel::BRep, analogous to OCCT TopoDS_Shape.
///
/// OCCT ref: TopoDS_Shape (TopoDS_Shape.hxx)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubShape {
    /// Reference to a face by its flat index in the BRep's face list.
    Face(usize),
    /// Reference to an edge by its index in the BRep's edge list.
    Edge(usize),
    /// Reference to a vertex by its index in the BRep's vertex list.
    Vertex(usize),
}

impl SubShape {
    /// Returns the entity type name as a string (e.g. "FACE", "EDGE", "VERTEX").
    pub fn type_name(&self) -> &'static str {
        match self {
            SubShape::Face(_) => "FACE",
            SubShape::Edge(_) => "EDGE",
            SubShape::Vertex(_) => "VERTEX",
        }
    }

    /// Returns the index of this shape in its entity array.
    pub fn index(&self) -> usize {
        match self {
            SubShape::Face(i) | SubShape::Edge(i) | SubShape::Vertex(i) => *i,
        }
    }
}

/// Statistics about the result of a boolean operation.
#[derive(Debug, Clone, Default)]
pub struct BooleanOpStatistics {
    /// Number of vertices in the result.
    pub n_vertices: usize,
    /// Number of edges in the result.
    pub n_edges: usize,
    /// Number of faces in the result.
    pub n_faces: usize,
    /// Number of shells in the result.
    pub n_shells: usize,
    /// Number of solids in the result.
    pub n_solids: usize,
}

/// `BuilderOperation` — a general boolean operation between two shapes.
///
/// This is the rcad equivalent of OCCT BRepAlgoAPI_BuilderOperation,
/// which is the base class for:
/// - BRepAlgoAPI_Common (intersection)
/// - BRepAlgoAPI_Fuse (union)
/// - BRepAlgoAPI_Cut (difference)
///
/// OCCT ref: BRepAlgoAPI_BuilderOperation (BRepAlgoAPI.cxx)
///
/// Wraps the rcad boolean pipeline (DS + PaveFiller + BooleanBuilder).
pub struct BuilderOperation {
    /// First input shape (object).
    shape_a: rcad_kernel::BRep,
    /// Second input shape (tool).
    shape_b: rcad_kernel::BRep,
    /// Boolean operation type.
    op_type: BooleanOpType,
    /// Tolerance for interference detection.
    tolerance: f64,
    /// Result shape (set after perform()).
    result: Option<rcad_kernel::BRep>,
    /// Error from the last perform() call.
    error: Option<BooleanError>,
}

impl BuilderOperation {
    /// Create a new BuilderOperation.
    ///
    /// OCCT ref: BRepAlgoAPI_Common/Fuse/Cut constructor (BRepAlgoAPI.cxx)
    ///
    /// The operation type determines which boolean is computed:
    /// - `Union` -> BRepAlgoAPI_Fuse
    /// - `Intersection` -> BRepAlgoAPI_Common
    /// - `Difference` -> BRepAlgoAPI_Cut
    pub fn new(a: rcad_kernel::BRep, b: rcad_kernel::BRep, op: BooleanOpType) -> Self {
        Self {
            shape_a: a,
            shape_b: b,
            op_type: op,
            tolerance: 1e-7,
            result: None,
            error: None,
        }
    }

    /// Set a custom tolerance for near-miss interference detection.
    /// Analogous to OCCT `BOPAlgo_Options::SetFuzzyValue()`.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol.max(1e-7);
    }

    /// Get the current tolerance value.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Perform the boolean operation.
    ///
    /// OCCT ref: BRepAlgoAPI_BuilderOperation::Build() (BRepAlgoAPI.cxx)
    ///
    /// Pipeline:
    /// 1. Build BOPDS_DS from the two shapes
    /// 2. Run BOPAlgo_PaveFiller to compute all intersections
    /// 3. Run BOPAlgo_Builder to build the result
    pub fn perform(&mut self) -> Result<&rcad_kernel::BRep, BooleanError> {
        self.result = None;
        self.error = None;

        // Check for empty inputs
        if self.shape_a.solids().is_empty() || self.shape_b.solids().is_empty() {
            let err = BooleanError::EmptyInput;
            self.error = Some(err);
            return Err(BooleanError::EmptyInput);
        }

        let a = &self.shape_a;
        let b = &self.shape_b;

        // Build BOPDS_DS (data structure) — OCCT: arguments are TopoDS_Solid/Shell directly
        let arg_a = Self::root_shape(a, 0);
        let arg_b = Self::root_shape(b, 1);
        let mut ds = DS::new();
        ds.set_arguments(vec![arg_a, arg_b]);
        ds.init(self.tolerance);

        // Run PaveFiller (BOPAlgo_PaveFiller::Perform)
        let mut filler = PaveFiller::new(&mut ds);
        filler.perform();
        let fuzz = filler.fuzzy_value();
        drop(filler);
        let mut builder = BooleanBuilder::new(&ds, self.op_type, fuzz);
        match builder.build() {
            Ok(brep) => {
                self.result = Some(brep);
                Ok(self.result.as_ref().unwrap())
            }
            Err(_) => {
                let err = BooleanError::InvalidOperation;
                self.error = Some(err);
                Err(BooleanError::InvalidOperation)
            }
        }
    }

    /// OCCT BOPAlgo_Algo::Prepare — extract root Solid/Shell from a BRep.
    fn root_shape(brep: &rcad_kernel::BRep, location: u32) -> Shape {
        for (i, ts) in brep.tshapes.iter().enumerate().rev() {
            match &**ts {
                TShape::Solid(_) | TShape::Shell(_) => {
                    return Shape::from_parts(ts.clone(), i, location, Orientation::Forward);
                }
                _ => {}
            }
        }
        // Fallback: should not happen for valid BReps
        Shape::new(Arc::new(TShape::Compound(
            brep.tshapes.iter().enumerate()
                .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, Orientation::Forward)).collect()
        )), location, Orientation::Forward)
    }

    /// Get the result shape.
    ///
    /// Panics if `perform()` has not been called or failed.
    pub fn shape(&self) -> &rcad_kernel::BRep {
        self.result
            .as_ref()
            .expect("perform() must be called before shape()")
    }

    /// Get the result shape, consuming the builder.
    pub fn into_shape(self) -> Option<rcad_kernel::BRep> {
        self.result
    }

    /// Get the error if the operation failed.
    pub fn error(&self) -> Option<&BooleanError> {
        self.error.as_ref()
    }

    /// Returns true if the operation has been performed successfully.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }

    // --- OCCT History Queries (stubs) ---
    // OCCT ref: BRepAlgoAPI_BuilderShape (BRepAlgoAPI_BuilderShape.hxx)

    /// Get all shapes that were modified during the operation.
    ///
    /// OCCT ref: BRepAlgoAPI_BuilderShape::Modified()
    pub fn modified(&self) -> Vec<SubShape> {
        // Stub: return nothing until history tracking is implemented
        Vec::new()
    }

    /// Get all shapes that were generated during the operation.
    ///
    /// OCCT ref: BRepAlgoAPI_BuilderShape::Generated()
    pub fn generated(&self) -> Vec<SubShape> {
        // Stub: return nothing until history tracking is implemented
        Vec::new()
    }

    /// Check if a source shape was deleted during the operation.
    ///
    /// OCCT ref: BRepAlgoAPI_BuilderShape::IsDeleted()
    pub fn is_deleted(&self, _source: &SubShape) -> bool {
        // Stub: return false until history tracking is implemented
        false
    }

    /// Get statistics about the result shape.
    pub fn statistics(&self) -> BooleanOpStatistics {
        let result = match self.result.as_ref() {
            Some(r) => r,
            None => return BooleanOpStatistics::default(),
        };
        let n_vertices = result.vertices().len();
        let n_edges = result.edges().len();
        let mut n_faces = 0;
        let mut n_shells = 0;
        let n_solids = result.solids().len();
        for solid in &result.solids() {
            for shell in &solid.shells {
                n_shells += 1;
                n_faces += shell.faces.len();
            }
        }
        BooleanOpStatistics {
            n_vertices,
            n_edges,
            n_faces,
            n_shells,
            n_solids,
        }
    }
}
