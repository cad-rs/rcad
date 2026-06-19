// OCCT-aligned: BOPAlgo_BuilderFace — face splitting using intersection curves.
//
// OCCT reference: BOPAlgo_BuilderFace.hxx, BOPAlgo_BuilderFace.cxx
//
// OCCT's BOPAlgo_BuilderFace is a separate class responsible for splitting a
// single face into sub-faces using intersection curves.  It inherits from
// BOPAlgo_BuilderArea and provides:
//
//   SetFace()        (BOPAlgo_BuilderFace.hxx L47)   — set the input face
//   Face()           (BOPAlgo_BuilderFace.hxx L50)   — get the input face
//   Perform()        (BOPAlgo_BuilderFace.hxx L53)   — main entry point
//     ├─ PerformShapesToAvoid()  (L62) — collect internal/same-orientation edges
//     ├─ PerformLoops()          (L69) — MakeBlocks + MakeWires
//     ├─ PerformAreas()          (L72) — classify wires as outer/holes
//     └─ PerformInternalShapes() (L75) — build finalized faces with internals
//   CheckData()      (L77) — validates input data
//
// In rcad this logic is embedded in BooleanBuilder methods and free functions
// in builder.rs.  BuilderFace provides an OCCT-aligned interface that delegates
// to the existing implementation.  Method bodies will be moved here incrementally.

use std::collections::HashMap;

use glam::DVec3;

use crate::bopds::ds::DS;
use crate::builder::{
    BooleanBuilder, WireFace, WireSegment, build_closed_wires, perform_areas,
};

/// OCCT-aligned: BOPAlgo_BuilderFace — splits a single face into sub-faces.
///
/// Wraps the face-splitting pipeline (edge collection, wire building, area
/// classification) behind an OCCT-aligned interface.  Currently delegates to
/// the existing pipeline in builder.rs.
///
/// OCCT reference: BOPAlgo_BuilderFace.hxx / BOPAlgo_BuilderFace.cxx
pub struct BuilderFace<'a> {
    /// The DS containing all topological entities
    ds: &'a DS,
    /// Index of the face to split
    face_idx: usize,
    /// Edge segments (boundary edges + intersection curves) — populated by perform()
    segments: Option<Vec<WireSegment>>,
    /// Resulting wire-faces (outer + inner wires) — populated by perform()
    wire_faces: Option<Vec<WireFace>>,
    /// Map from canonical vertex index to 3D position — populated by perform()
    vertex_positions: Option<HashMap<usize, DVec3>>,
}

impl<'a> BuilderFace<'a> {
    /// OCCT-aligned: SetFace equivalent via constructor.
    ///
    /// Stores the DS reference and the face index to operate on.
    /// In OCCT this is SetFace(const TopoDS_Face&) (BOPAlgo_BuilderFace.hxx L47).
    pub fn new(ds: &'a DS, face_idx: usize) -> Self {
        Self {
            ds,
            face_idx,
            segments: None,
            wire_faces: None,
            vertex_positions: None,
        }
    }

    /// OCCT-aligned: Perform — main entry point for face splitting.
    ///
    /// Delegates to BooleanBuilder::split_face_occt_wire_pipeline, which
    /// implements the OCCT flow:
    ///
    /// 1. Collect face edge segments (boundary + intersection curves)
    /// 2. Build closed wires from blocks (MakeBlocks + MakeWires)
    /// 3. Classify wires as outer/holes (PerformAreas)
    ///
    /// OCCT reference: BOPAlgo_BuilderFace.hxx L53, BOPAlgo_BuilderFace.cxx L58-149
    ///
    /// Returns true on success, false if the face produced no sub-faces.
    ///
    /// ✅ OCCT-aligned: delegates to BooleanBuilder which follows OCCT's
    ///    Perform → PerformanceLoops → PerformAreas flow.
    /// ⏳ Partial alignment: the internal logic still lives in BooleanBuilder
    ///    rather than in this struct.
    pub fn perform(&mut self, builder: &BooleanBuilder) -> bool {
        if let Some((segments, wfs, vertex_positions)) =
            builder.split_face_occt_wire_pipeline(self.face_idx)
        {
            self.segments = Some(segments);
            self.wire_faces = Some(wfs);
            self.vertex_positions = Some(vertex_positions);
            true
        } else {
            false
        }
    }

    /// OCCT-aligned: PerformLoops / MakeWires — build closed wires from edge blocks.
    ///
    /// Calls the free function `build_closed_wires()` which implements:
    /// - MakeConnexityBlocks: BFS grouping by shared vertices
    /// - Regular blocks (degree=2): simple walk
    /// - Irregular blocks (degree>2): SmartMap + angle-based Path walking
    ///
    /// OCCT reference: BOPAlgo_BuilderFace.hxx L69 (PerformLoops),
    ///                 BOPAlgo_BuilderFace.cxx L239-606
    ///
    /// ✅ OCCT-aligned: delegates to build_closed_wires which follows the
    ///    OCCT MakeBlocks + MakeWires flow.
    pub fn build_wires(
        &self,
        segments: &mut Vec<WireSegment>,
    ) -> (
        Vec<Vec<usize>>,
        Vec<Vec<usize>>,
        HashMap<usize, DVec3>,
    ) {
        build_closed_wires(segments, self.ds, self.face_idx)
    }

    /// OCCT-aligned: PerformAreas — classify wires as outer boundary or hole.
    ///
    /// Calls the free function `perform_areas()` which:
    /// - Computes 3D boundary polygon and centroid for each wire
    /// - Sorts by projected area (largest = potential outer)
    /// - Classifies wires as growth (outer) or hole via point-in-polygon
    /// - Assigns holes to enclosing growths
    ///
    /// OCCT reference: BOPAlgo_BuilderFace.hxx L72 (PerformAreas),
    ///                 BOPAlgo_BuilderFace.cxx L428-613
    ///
    /// ✅ OCCT-aligned: delegates to perform_areas which follows the
    ///    OCCT area classification flow.
    pub fn perform_areas(
        &self,
        wires: &[Vec<usize>],
        internal_wires: &[Vec<usize>],
        segments: &[WireSegment],
    ) -> Vec<WireFace> {
        perform_areas(wires, internal_wires, segments, self.ds, self.face_idx)
    }

    // ── Accessors ──────────────────────────────────────────────────────────

    /// Returns the face index this builder operates on.
    pub fn face_idx(&self) -> usize {
        self.face_idx
    }

    /// Returns the edge segments after perform() completes, if available.
    pub fn segments(&self) -> Option<&[WireSegment]> {
        self.segments.as_deref()
    }

    /// Returns the wire-faces after perform() completes, if available.
    pub fn wire_faces(&self) -> Option<&[WireFace]> {
        self.wire_faces.as_deref()
    }

    /// Returns the vertex position map after perform() completes, if available.
    pub fn vertex_positions(&self) -> Option<&HashMap<usize, DVec3>> {
        self.vertex_positions.as_ref()
    }

    /// Consume the builder and return the wire-faces, if any.
    pub fn into_wire_faces(self) -> Option<Vec<WireFace>> {
        self.wire_faces
    }

    /// Consume the builder and return the segments, if any.
    pub fn into_segments(self) -> Option<Vec<WireSegment>> {
        self.segments
    }

    /// Consume the builder and return the vertex positions, if any.
    pub fn into_vertex_positions(self) -> Option<HashMap<usize, DVec3>> {
        self.vertex_positions
    }
}
