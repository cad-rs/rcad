//! OCCT BRepAlgoAPI_BuilderOperation equivalent 鈥?boolean operation wrapper with history.
//!
//! Provides a general-purpose boolean operation struct analogous to
//! OCCT BRepAlgoAPI_BuilderOperation (base class for BRepAlgoAPI_Common/Fuse/Cut).
//!
//! OCCT references:
//! - BRepAlgoAPI_Common / BRepAlgoAPI_Fuse / BRepAlgoAPI_Cut (BRepAlgoAPI.cxx)
//! - BRepAlgoAPI_BuilderOperation 鈫?BRepAlgoAPI_BuilderShape 鈫?BRepAlgoAPI_Algo
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::brep_algo_api::builder_operation::BooleanOp;
//! use rcad_algorithms::BooleanOpType;
//! use rcad_kernel::{BRep, PrimitiveSolid};
//!
//! let a = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
//! let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 3.0 });
//!
//! let mut op = BooleanOp::new(a, b, BooleanOpType::Union);
//! if let Ok(result) = op.perform() {
//!     println!("Result has {} faces", op.statistics().n_faces);
//! }
//! ```


use crate::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::bopds::ds::DS;
use crate::bvh::Bvh;
use crate::history::BooleanHistory;
use crate::pave_filler::PaveFiller;
use crate::tolerance::TOLERANCE_ABS;
use rcad_kernel::topods;

/// A reference to a sub-shape in a rcad_kernel::BRep, analogous to OCCT TopoDS_Shape.
///
/// OCCT ref: TopoDS_Shape (TopoDS_Shape.hxx)
///
/// This is a lightweight handle that identifies a specific sub-shape
/// (face, edge, or vertex) in a BRep by its index in the corresponding array.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShapeRef {
    /// Reference to a face by its flat index in the BRep's face list.
    Face(usize),
    /// Reference to an edge by its index in the BRep's edge list.
    Edge(usize),
    /// Reference to a vertex by its index in the BRep's vertex list.
    Vertex(usize),
}

impl ShapeRef {
    /// Returns the entity type name as a string (e.g. "FACE", "EDGE", "VERTEX").
    ///
    /// OCCT ref: TopoDS_Shape::ShapeType() string mapping.
    pub fn type_name(&self) -> &'static str {
        match self {
            ShapeRef::Face(_) => "FACE",
            ShapeRef::Edge(_) => "EDGE",
            ShapeRef::Vertex(_) => "VERTEX",
        }
    }

    /// Returns the index of this shape in its entity array.
    pub fn index(&self) -> usize {
        match self {
            ShapeRef::Face(i) | ShapeRef::Edge(i) | ShapeRef::Vertex(i) => *i,
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

/// BooleanOp 鈥?a general boolean operation between two shapes with history tracking.
///
/// This is the rcad equivalent of OCCT BRepAlgoAPI_BuilderOperation,
/// which is the base class for:
/// - BRepAlgoAPI_Common (intersection)
/// - BRepAlgoAPI_Fuse (union)
/// - BRepAlgoAPI_Cut (difference)
///
/// OCCT ref: BRepAlgoAPI_BuilderOperation (BRepAlgoAPI.cxx)
///
/// Wraps the rcad boolean pipeline (DS + PaveFiller + BooleanBuilder) in
/// a struct that tracks shape history (Modified/Generated/IsDeleted).
pub struct BooleanOp {
    /// First input shape (object).
    shape_a: rcad_kernel::BRep,
    /// Second input shape (tool).
    shape_b: rcad_kernel::BRep,
    /// Boolean operation type.
    op_type: BooleanOpType,
    /// History of the operation (set after perform()).
    history: Option<BooleanHistory>,
    /// Tolerance for interference detection.
    tolerance: f64,
    /// Result shape (set after perform()).
    result: Option<rcad_kernel::BRep>,
    /// Error from the last perform() call.
    error: Option<BooleanError>,
}

impl BooleanOp {
    /// Create a new BooleanOp.
    ///
    /// OCCT ref: BRepAlgoAPI_Common/Fuse/Cut constructor (BRepAlgoAPI.cxx)
    ///
    /// The operation type determines which boolean is computed:
    /// - `Union` 鈫?BRepAlgoAPI_Fuse
    /// - `Intersection` 鈫?BRepAlgoAPI_Common
    /// - `Difference` 鈫?BRepAlgoAPI_Cut
    pub fn new(a: rcad_kernel::BRep, b: rcad_kernel::BRep, op: BooleanOpType) -> Self {
        Self {
            shape_a: a,
            shape_b: b,
            op_type: op,
            history: None,
            tolerance: TOLERANCE_ABS,
            result: None,
            error: None,
        }
    }

    /// Set a custom tolerance for near-miss interference detection.
    ///
    /// Must be >= TOLERANCE_ABS. Values below are clamped.
    /// Analogous to OCCT `BOPAlgo_Options::SetFuzzyValue()`.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol.max(TOLERANCE_ABS);
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
    /// 1. 鉁?Build BOPDS_DS from the two shapes
    /// 2. 鉁?Run BOPAlgo_PaveFiller to compute all intersections
    /// 3. 鉁?Build container images (FillImagesContainers)
    /// 4. 鉁?Run BOPAlgo_Builder to build the result
    ///
    /// Returns a reference to the result BRep on success.
    pub fn perform(&mut self) -> Result<&rcad_kernel::BRep, BooleanError> {
        self.result = None;
        self.error = None;
        self.history = None;

        // 鉁?OCCT-aligned: BOPAlgo_BOP::CheckInputData 鈥?verify valid inputs
        if self.shape_a.solids.is_empty() || self.shape_b.solids.is_empty() {
            let err = BooleanError::EmptyInput;
            self.error = Some(err);
            return Err(BooleanError::EmptyInput);
        }

        // Ensure geometry is populated
        let a = self.ensure_geometry(&self.shape_a);
        let b = self.ensure_geometry(&self.shape_b);

        // 鉁?OCCT-aligned: Build BOPDS_DS (data structure) from the two shapes
        // OCCT ref: BOPAlgo_PaveFiller::Perform 鈫?BOPDS_DS::Alloc
        let a_t = a.to_topods();
        let b_t = b.to_topods();
        let mut ds = DS::new_from_topods(&a_t, &b_t, self.tolerance.max(TOLERANCE_ABS));

        // 鉁?OCCT-aligned: Build BVH acceleration (optional in OCCT)
        let bvh_a = Bvh::build(&a);
        let bvh_b = Bvh::build(&b);

        // 鉁?OCCT-aligned: Run PaveFiller (BOPAlgo_PaveFiller::Perform)
        let mut brep = rcad_kernel::topods::BRep::new();
        let (face_refs, ic_edge_map) = {
            let mut filler = PaveFiller::with_bvh_and_brep(&mut ds, &bvh_a, &bvh_b, &mut brep);
            filler.perform();
            (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
        };

        // 鉁?OCCT-aligned: FillImagesContainers 鈥?build wire/shell images
        ds.build_container_images();

        // 鉁?OCCT-aligned: Build result
        let builder = BooleanBuilder::with_brep(&ds, self.op_type, brep, face_refs, ic_edge_map);
        let (t, bool_history) = builder.build_with_history()?;

        self.history = Some(bool_history);
        self.result = Some(rcad_kernel::BRep::from_topods(&t));
        Ok(self.result.as_ref().unwrap())
    }

    /// Ensure geometry is populated for primitive shapes.
    fn ensure_geometry(&self, brep: &rcad_kernel::BRep) -> rcad_kernel::BRep {
        if brep.geom.surfaces.is_empty() && !brep.solids.is_empty() {
            let mut result = brep.clone();
            crate::geom_populate::populate_box_geom(&mut result);
            result
        } else {
            brep.clone()
        }
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

    // 鈹€鈹€ OCCT History Queries 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    //
    // OCCT ref: BRepAlgoAPI_BuilderShape (BRepAlgoAPI_BuilderShape.hxx)
    //
    // These methods provide OCCT-compatible shape history tracking.
    // 鉁?OCCT-aligned: stubs return all known history.
    // 鈴?Side-distinction (A vs B) not implemented 鈥?section-style operations
    //     always involve both shapes; HasAncestorFaceOn1/2 always true.

    /// Get all shapes that were modified (split or changed) during the operation.
    ///
    /// OCCT ref: BRepAlgoAPI_BuilderShape::Modified()
    ///
    /// Returns all result shapes whose origin traces back to an input shape
    /// (i.e., faces/edges/vertices that were split or carried through).
    ///
    /// 鉁?OCCT-aligned: returns all modified shapes.
    pub fn modified(&self) -> Vec<ShapeRef> {
        let Some(ref h) = self.history else {
            return Vec::new();
        };
        let mut result = Vec::new();

        // Collect modified faces
        for (res_idx, origin) in h.face_origins.iter().enumerate() {
            match origin {
                crate::history::FaceOrigin::FromA(_) | crate::history::FaceOrigin::FromB(_) => {
                    result.push(ShapeRef::Face(res_idx));
                }
                _ => {}
            }
        }
        for &(res_idx, ref origin) in &h.co_face_origins {
            match origin {
                crate::history::FaceOrigin::FromA(_) | crate::history::FaceOrigin::FromB(_) => {
                    result.push(ShapeRef::Face(res_idx));
                }
                _ => {}
            }
        }

        // Collect modified edges
        for (res_idx, origin) in h.edge_origins.iter().enumerate() {
            match origin {
                crate::history::EdgeOrigin::FromA(_)
                | crate::history::EdgeOrigin::FromB(_)
                | crate::history::EdgeOrigin::SplitFromA(_)
                | crate::history::EdgeOrigin::SplitFromB(_) => {
                    result.push(ShapeRef::Edge(res_idx));
                }
                _ => {}
            }
        }

        // Collect modified vertices
        for (res_idx, origin) in h.vertex_origins.iter().enumerate() {
            match origin {
                crate::history::VertexOrigin::FromA(_)
                | crate::history::VertexOrigin::FromB(_) => {
                    result.push(ShapeRef::Vertex(res_idx));
                }
                _ => {}
            }
        }

        result
    }

    /// Get all shapes that were generated (newly created) during the operation.
    ///
    /// OCCT ref: BRepAlgoAPI_BuilderShape::Generated()
    ///
    /// Returns result shapes that were created as a direct result of the
    /// boolean operation (intersection edges, vertices, etc.).
    ///
    /// 鉁?OCCT-aligned: returns generated faces/edges/vertices.
    pub fn generated(&self) -> Vec<ShapeRef> {
        let Some(ref h) = self.history else {
            return Vec::new();
        };
        let mut result = Vec::new();

        // Collect generated faces
        for (res_idx, origin) in h.face_origins.iter().enumerate() {
            if matches!(origin, crate::history::FaceOrigin::Generated) {
                result.push(ShapeRef::Face(res_idx));
            }
        }
        for &(res_idx, ref origin) in &h.co_face_origins {
            if matches!(origin, crate::history::FaceOrigin::Generated) {
                result.push(ShapeRef::Face(res_idx));
            }
        }

        // Collect generated edges
        for (res_idx, origin) in h.edge_origins.iter().enumerate() {
            if matches!(origin, crate::history::EdgeOrigin::Generated) {
                result.push(ShapeRef::Edge(res_idx));
            }
        }

        // Collect intersection vertices
        for (res_idx, origin) in h.vertex_origins.iter().enumerate() {
            if matches!(origin, crate::history::VertexOrigin::Intersection) {
                result.push(ShapeRef::Vertex(res_idx));
            }
        }

        result
    }

    /// Check if a source shape was deleted during the operation.
    ///
    /// OCCT ref: BRepAlgoAPI_BuilderShape::IsDeleted(const TopoDS_Shape&)
    ///
    /// Returns true if the source shape was entirely removed.
    /// For faces: checks that the source face index appears in the deleted list.
    /// For edges/vertices: checks the deletion tracker.
    ///
    /// 鉁?OCCT-aligned: face deletion tracking via history.
    pub fn is_deleted(&self, source: &ShapeRef) -> bool {
        let Some(ref h) = self.history else {
            return false;
        };
        match source {
            ShapeRef::Face(idx) => {
                h.deleted_from_a.contains(idx) || h.deleted_from_b.contains(idx)
            }
            ShapeRef::Edge(idx) => {
                h.tracker.deleted_edges().any(|e| e.entity_index == *idx)
            }
            ShapeRef::Vertex(idx) => {
                h.tracker.deleted_vertices().any(|v| v.entity_index == *idx)
            }
        }
    }

    /// Get statistics about the result shape.
    ///
    /// Returns counts of vertices, edges, faces, shells, and solids
    /// in the result BRep. Returns default (all zeros) if the operation
    /// has not been performed yet.
    pub fn statistics(&self) -> BooleanOpStatistics {
        let result = match self.result.as_ref() {
            Some(r) => r,
            None => return BooleanOpStatistics::default(),
        };
        let n_vertices = result.vertices.len();
        let n_edges = result.edges.len();
        let mut n_faces = 0;
        let mut n_shells = 0;
        let n_solids = result.solids.len();
        for solid in &result.solids {
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


