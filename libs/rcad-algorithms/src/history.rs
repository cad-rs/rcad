use rcad_kernel::{BRep, PersistentNamingHooks, TopoEntityRef};

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

/// Report produced when propagating persistent names through a boolean result.
///
/// This captures which source names could not be mapped into the result and
/// which names had to be duplicated because a single source entity generated
/// multiple result entities.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BooleanNamingPropagationReport {
    pub dropped_from_a: Vec<String>,
    pub dropped_from_b: Vec<String>,
    pub duplicated_from_a: Vec<String>,
    pub duplicated_from_b: Vec<String>,
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

    /// Propagate source persistent names through this boolean history into the
    /// `result_brep` topology.
    ///
    /// Face, edge, vertex, and baseline solid names from both inputs are mapped
    /// onto result entities according to `face_origins`, `edge_origins`,
    /// `vertex_origins`, and `solid_origins`.
    ///
    /// When one source entity produces multiple result entities (for example a
    /// split edge), the original name is bound to the first result entity and
    /// deterministic suffixed variants (`name@1`, `name@2`, ...) are bound to the
    /// remaining ones.
    pub fn propagate_persistent_naming(
        &self,
        result_brep: &BRep,
        names_a: &PersistentNamingHooks,
        names_b: &PersistentNamingHooks,
    ) -> (PersistentNamingHooks, BooleanNamingPropagationReport) {
        let mut out = PersistentNamingHooks::new();
        let mut report = BooleanNamingPropagationReport::default();

        self.propagate_from_source(
            &mut out,
            names_a,
            InputSide::A,
            &mut report.dropped_from_a,
            &mut report.duplicated_from_a,
        );
        self.propagate_from_source(
            &mut out,
            names_b,
            InputSide::B,
            &mut report.dropped_from_b,
            &mut report.duplicated_from_b,
        );

        out.retain_valid_for_brep(result_brep);
        (out, report)
    }

    fn propagate_from_source(
        &self,
        out: &mut PersistentNamingHooks,
        source: &PersistentNamingHooks,
        side: InputSide,
        dropped: &mut Vec<String>,
        duplicated: &mut Vec<String>,
    ) {
        for (name, target_ref) in source.iter() {
            let matches = self.matching_result_entities(side, target_ref);
            if matches.is_empty() {
                dropped.push(name.to_string());
                continue;
            }
            if matches.len() > 1 {
                duplicated.push(name.to_string());
            }
            for (suffix_idx, result_ref) in matches.into_iter().enumerate() {
                let bound_name = if suffix_idx == 0 {
                    name.to_string()
                } else {
                    format!("{name}@{suffix_idx}")
                };
                bind_unique(out, bound_name, result_ref);
            }
        }
    }

    fn matching_result_entities(&self, side: InputSide, target_ref: TopoEntityRef) -> Vec<TopoEntityRef> {
        match target_ref {
            TopoEntityRef::Face(source_idx) => self
                .face_origins
                .iter()
                .enumerate()
                .filter_map(|(result_idx, origin)| match (side, origin) {
                    (InputSide::A, FaceOrigin::FromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Face(result_idx))
                    }
                    (InputSide::B, FaceOrigin::FromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Face(result_idx))
                    }
                    _ => None,
                })
                .collect(),
            TopoEntityRef::Edge(source_idx) => self
                .edge_origins
                .iter()
                .enumerate()
                .filter_map(|(result_idx, origin)| match (side, origin) {
                    (InputSide::A, EdgeOrigin::FromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    (InputSide::A, EdgeOrigin::SplitFromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    (InputSide::B, EdgeOrigin::FromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    (InputSide::B, EdgeOrigin::SplitFromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Edge(result_idx))
                    }
                    _ => None,
                })
                .collect(),
            TopoEntityRef::Vertex(source_idx) => self
                .vertex_origins
                .iter()
                .enumerate()
                .filter_map(|(result_idx, origin)| match (side, origin) {
                    (InputSide::A, VertexOrigin::FromA(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Vertex(result_idx))
                    }
                    (InputSide::B, VertexOrigin::FromB(idx)) if *idx == source_idx => {
                        Some(TopoEntityRef::Vertex(result_idx))
                    }
                    _ => None,
                })
                .collect(),
            TopoEntityRef::Solid(source_idx) => {
                if source_idx != 0 {
                    return Vec::new();
                }
                self.solid_origins
                    .iter()
                    .enumerate()
                    .filter_map(|(result_idx, origin)| match (side, origin) {
                        (InputSide::A, SolidOrigin::FromA) => Some(TopoEntityRef::Solid(result_idx)),
                        (InputSide::B, SolidOrigin::FromB) => Some(TopoEntityRef::Solid(result_idx)),
                        _ => None,
                    })
                    .collect()
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputSide {
    A,
    B,
}

fn bind_unique(out: &mut PersistentNamingHooks, preferred_name: String, target_ref: TopoEntityRef) {
    if out.resolve(&preferred_name).is_none() {
        out.bind(preferred_name, target_ref);
        return;
    }
    let mut suffix_idx = 1usize;
    loop {
        let candidate = format!("{preferred_name}@{suffix_idx}");
        if out.resolve(&candidate).is_none() {
            out.bind(candidate, target_ref);
            return;
        }
        suffix_idx += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{PrimitiveSolid, TopoEntityRef};

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    #[test]
    fn propagate_persistent_naming_maps_face_edge_vertex_and_solid_origins() {
        let result_brep = unit_box();
        let history = BooleanHistory {
            face_origins: vec![FaceOrigin::FromA(0), FaceOrigin::FromB(1)],
            edge_origins: vec![EdgeOrigin::FromA(0), EdgeOrigin::SplitFromA(0), EdgeOrigin::FromB(1)],
            vertex_origins: vec![VertexOrigin::FromA(0), VertexOrigin::Intersection, VertexOrigin::FromB(1)],
            shell_origins: vec![],
            solid_origins: vec![SolidOrigin::FromA],
        };

        let mut names_a = PersistentNamingHooks::new();
        names_a.bind("face_a", TopoEntityRef::Face(0));
        names_a.bind("edge_a", TopoEntityRef::Edge(0));
        names_a.bind("vertex_a", TopoEntityRef::Vertex(0));
        names_a.bind("solid_a", TopoEntityRef::Solid(0));

        let mut names_b = PersistentNamingHooks::new();
        names_b.bind("face_b", TopoEntityRef::Face(1));
        names_b.bind("edge_b", TopoEntityRef::Edge(1));
        names_b.bind("vertex_b", TopoEntityRef::Vertex(1));
        names_b.bind("solid_b", TopoEntityRef::Solid(0));

        let (result_names, report) = history.propagate_persistent_naming(&result_brep, &names_a, &names_b);

        assert_eq!(result_names.resolve("face_a"), Some(TopoEntityRef::Face(0)));
        assert_eq!(result_names.resolve("face_b"), Some(TopoEntityRef::Face(1)));
        assert_eq!(result_names.resolve("edge_a"), Some(TopoEntityRef::Edge(0)));
        assert_eq!(result_names.resolve("edge_a@1"), Some(TopoEntityRef::Edge(1)));
        assert_eq!(result_names.resolve("edge_b"), Some(TopoEntityRef::Edge(2)));
        assert_eq!(result_names.resolve("vertex_a"), Some(TopoEntityRef::Vertex(0)));
        assert_eq!(result_names.resolve("vertex_b"), Some(TopoEntityRef::Vertex(2)));
        assert_eq!(result_names.resolve("solid_a"), Some(TopoEntityRef::Solid(0)));
        assert_eq!(result_names.resolve("solid_b"), None);

        assert!(report.dropped_from_a.is_empty());
        assert_eq!(report.dropped_from_b, vec!["solid_b".to_string()]);
        assert_eq!(report.duplicated_from_a, vec!["edge_a".to_string()]);
        assert!(report.duplicated_from_b.is_empty());
    }

    #[test]
    fn propagate_persistent_naming_disambiguates_cross_input_name_collisions() {
        let result_brep = unit_box();
        let history = BooleanHistory {
            face_origins: vec![FaceOrigin::FromA(0), FaceOrigin::FromB(0)],
            edge_origins: vec![],
            vertex_origins: vec![],
            shell_origins: vec![],
            solid_origins: vec![],
        };

        let mut names_a = PersistentNamingHooks::new();
        names_a.bind("shared_face", TopoEntityRef::Face(0));
        let mut names_b = PersistentNamingHooks::new();
        names_b.bind("shared_face", TopoEntityRef::Face(0));

        let (result_names, report) = history.propagate_persistent_naming(&result_brep, &names_a, &names_b);

        assert_eq!(result_names.resolve("shared_face"), Some(TopoEntityRef::Face(0)));
        assert_eq!(result_names.resolve("shared_face@1"), Some(TopoEntityRef::Face(1)));
        assert!(report.dropped_from_a.is_empty());
        assert!(report.dropped_from_b.is_empty());
    }
}
