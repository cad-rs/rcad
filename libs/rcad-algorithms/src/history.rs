/// Tracks the origin of each face in a boolean operation result.
///
/// Analogous to OCCT `BRepAlgoAPI_BuilderShape::Modified/Generated/Deleted`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceOrigin {
    /// Face came from solid A; value is the DS face index of the source face.
    FromA(usize),
    /// Face came from solid B; value is the DS face index of the source face.
    FromB(usize),
    /// Face was generated at the intersection boundary (not yet produced by this
    /// implementation — reserved for future use).
    Generated,
}

/// Per-face origin map for a boolean operation result.
///
/// `face_origins[i]` gives the origin of `result_brep.solids[0].shells[0].faces[i]`.
#[derive(Debug, Clone)]
pub struct BooleanHistory {
    pub face_origins: Vec<FaceOrigin>,
}

impl BooleanHistory {
    /// Returns the origin of face `idx` in the result BRep.
    pub fn face_origin(&self, idx: usize) -> FaceOrigin {
        self.face_origins[idx]
    }

    /// Number of result faces tracked.
    pub fn len(&self) -> usize {
        self.face_origins.len()
    }

    pub fn is_empty(&self) -> bool {
        self.face_origins.is_empty()
    }

    /// How many result faces came from solid A.
    pub fn count_from_a(&self) -> usize {
        self.face_origins
            .iter()
            .filter(|o| matches!(o, FaceOrigin::FromA(_)))
            .count()
    }

    /// How many result faces came from solid B.
    pub fn count_from_b(&self) -> usize {
        self.face_origins
            .iter()
            .filter(|o| matches!(o, FaceOrigin::FromB(_)))
            .count()
    }
}
