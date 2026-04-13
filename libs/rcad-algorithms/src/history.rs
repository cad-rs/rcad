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

/// Tracks the origin of each edge in a boolean operation result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeOrigin {
    /// Edge came from solid A; value is the original edge index in solid A.
    FromA(usize),
    /// Edge came from solid B; value is the original edge index in solid B.
    FromB(usize),
    /// Edge is an intersection edge generated at the boolean boundary.
    Generated,
    /// Edge was created from a partial (split) segment of an original edge in A.
    SplitFromA(usize),
    /// Edge was created from a partial (split) segment of an original edge in B.
    SplitFromB(usize),
}

/// Tracks the origin of each vertex in a boolean operation result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexOrigin {
    /// Vertex came from solid A; value is the original vertex index in solid A.
    FromA(usize),
    /// Vertex came from solid B; value is the original vertex index in solid B.
    FromB(usize),
    /// Vertex was created at an A-B intersection point.
    Intersection,
}

/// Aggregate origin of a result shell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellOrigin {
    /// Every tracked face in the shell came from solid A.
    FromA,
    /// Every tracked face in the shell came from solid B.
    FromB,
    /// Every tracked face in the shell was generated.
    Generated,
    /// The shell contains a mixture of A/B/generated faces.
    Mixed,
}

/// Aggregate origin of a result solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidOrigin {
    /// Every tracked shell in the solid came from solid A.
    FromA,
    /// Every tracked shell in the solid came from solid B.
    FromB,
    /// Every tracked shell in the solid was generated.
    Generated,
    /// The solid contains a mixture of A/B/generated shells or faces.
    Mixed,
}

/// Per-face origin map for a boolean operation result.
///
/// `face_origins[i]` gives the origin of `result_brep.solids[0].shells[0].faces[i]`.
#[derive(Debug, Clone)]
pub struct BooleanHistory {
    pub face_origins: Vec<FaceOrigin>,
    /// Per-edge origin map. `edge_origins[i]` gives the origin of `result_brep.edges[i]`.
    ///
    /// Empty when edge history was not requested (standard `boolean_op_with_history` path
    /// does not yet populate this; use `boolean_op_with_full_history` for edge tracking).
    pub edge_origins: Vec<EdgeOrigin>,
    /// Per-vertex origin map. `vertex_origins[i]` gives the origin of `result_brep.vertices[i]`.
    ///
    /// Empty when vertex history was not requested.
    pub vertex_origins: Vec<VertexOrigin>,
    /// Per-shell aggregate origin map. Flattened in the same order as `result_brep.solids[*].shells[*]`.
    pub shell_origins: Vec<ShellOrigin>,
    /// Per-solid aggregate origin map. `solid_origins[i]` gives the origin of `result_brep.solids[i]`.
    pub solid_origins: Vec<SolidOrigin>,
}

impl BooleanHistory {
    /// Returns the origin of face `idx` in the result BRep.
    pub fn face_origin(&self, idx: usize) -> FaceOrigin {
        self.face_origins[idx]
    }

    /// Returns the origin of edge `idx` in the result BRep.
    /// Returns `None` if edge history was not recorded.
    pub fn edge_origin(&self, idx: usize) -> Option<EdgeOrigin> {
        self.edge_origins.get(idx).copied()
    }

    /// Returns the origin of vertex `idx` in the result BRep.
    /// Returns `None` if vertex history was not recorded.
    pub fn vertex_origin(&self, idx: usize) -> Option<VertexOrigin> {
        self.vertex_origins.get(idx).copied()
    }

    /// Returns the aggregate origin of shell `idx` in the flattened result BRep.
    pub fn shell_origin(&self, idx: usize) -> Option<ShellOrigin> {
        self.shell_origins.get(idx).copied()
    }

    /// Returns the aggregate origin of solid `idx` in the result BRep.
    pub fn solid_origin(&self, idx: usize) -> Option<SolidOrigin> {
        self.solid_origins.get(idx).copied()
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

    /// How many result edges came from solid A (including splits).
    pub fn edge_count_from_a(&self) -> usize {
        self.edge_origins
            .iter()
            .filter(|o| matches!(o, EdgeOrigin::FromA(_) | EdgeOrigin::SplitFromA(_)))
            .count()
    }

    /// How many result edges came from solid B (including splits).
    pub fn edge_count_from_b(&self) -> usize {
        self.edge_origins
            .iter()
            .filter(|o| matches!(o, EdgeOrigin::FromB(_) | EdgeOrigin::SplitFromB(_)))
            .count()
    }

    /// How many result edges were generated at the intersection.
    pub fn edge_count_generated(&self) -> usize {
        self.edge_origins
            .iter()
            .filter(|o| matches!(o, EdgeOrigin::Generated))
            .count()
    }

    /// How many result solids contain contributions from both inputs and/or generated topology.
    pub fn solid_count_mixed(&self) -> usize {
        self.solid_origins
            .iter()
            .filter(|o| matches!(o, SolidOrigin::Mixed))
            .count()
    }
}
