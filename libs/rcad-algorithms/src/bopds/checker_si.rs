/// OCCT-aligned CheckerSI: validates a single shape for self-interference.
///
/// OCCT reference: BOPAlgo_CheckerSI (BOPAlgo_CheckerSI.cxx L1-300).
///
/// Reuses the existing PaveFiller to detect self-interferences within
/// a single shape by loading it as both inputs and running the standard
/// interference detection pipeline, then filtering out trivial A-B pairs
/// where the entities are the same original entity.
///
/// # Usage
///
/// ```ignore
/// let mut checker = CheckerSI::new();
/// checker.set_level_of_check(3);
/// checker.perform(&brep);
/// if checker.has_interferences() {
///     for interf in checker.get_interferences() {
///         // handle self-interference
///     }
/// }
/// ```

use crate::bopds::ds::{DS, Interference};
use crate::pave_filler::PaveFiller;
use rcad_kernel::topods;

/// CheckerSI validates a single shape for self-interference.
///
/// OCCT ref: BOPAlgo_CheckerSI (BOPAlgo_CheckerSI.cxx).
///
/// Reuses PaveFiller with the shape loaded as both inputs.
/// After `Perform()`, filtered interferences represent true
/// self-interferences within the original shape.
pub struct CheckerSI {
    /// Level of check (0-9), controlling which interference types are reported.
    ///
    /// OCCT BOPAlgo_CheckerSI::SetLevelOfCheck:
    /// - Level 0: VE, EE, EF, FF (basic entity checks)
    /// - Level 1: add VV
    /// - Level 2: add VF
    /// - Level 3+: all
    level_of_check: u8,
    /// Interferences detected in the last Perform() call, filtered by
    /// non-triviality and level of check.
    interferences: Vec<Interference>,
    /// Whether self-interferences were found in the last Perform().
    has_interferences: bool,
}

impl CheckerSI {
    /// Create a new CheckerSI with default LevelOfCheck = 0.
    pub fn new() -> Self {
        Self {
            level_of_check: 0,
            interferences: Vec::new(),
            has_interferences: false,
        }
    }

    /// Set the level of check (clamped to 0-9).
    ///
    /// OCCT ref: BOPAlgo_CheckerSI::SetLevelOfCheck (BOPAlgo_CheckerSI.cxx L50-70).
    ///
    /// Controls which intersection passes are included:
    /// - Level 0: VE, EE, EF, FF
    /// - Level 1: +VV
    /// - Level 2: +VF
    /// - Level 3+: all
    pub fn set_level_of_check(&mut self, level: u8) {
        self.level_of_check = level.min(9);
    }

    /// Get the current level of check.
    pub fn level_of_check(&self) -> u8 {
        self.level_of_check
    }

    /// Run the self-interference check on a single topods::BRep.
    ///
    /// OCCT ref: BOPAlgo_CheckerSI::Perform (BOPAlgo_CheckerSI.cxx L80-200).
    ///
    /// Loads the BRep as both ShapeA and ShapeB in the DS, then runs PaveFiller.
    /// After filtering trivial same-entity pairs and applying LevelOfCheck,
    /// the remaining interferences represent true self-interferences.
    pub fn perform(&mut self, brep: &topods::BRep) {
        // Load the single shape twice (as both A and B) into the DS.
        // PaveFiller compares A entities against B entities, so loading the
        // same shape twice makes it detect intersections between any pair
        // of distinct sub-shapes within the original shape.
        let mut ds = DS::new_empty();
        // Run the PaveFiller interference detection pipeline.
        // This executes all six intersection passes (VV/VE/EE/VF/EF/FF)
        // plus edge splitting and common block detection.
        let mut pf = PaveFiller::new(&mut ds);
        pf.perform(brep, brep);
        // a_vertex/edge/face counts set by init() inside perform()
        let a_vc = ds.a_vertex_count;
        let a_ec = ds.a_edge_count;
        let a_fc = ds.a_face_count;

        // Filter interferences:
        // 1. Remove trivial A-B pairs (same original entity from two copies).
        // 2. Apply LevelOfCheck filter.
        let level = self.level_of_check;
        let mut filtered: Vec<Interference> = Vec::new();
        for inf in &ds.interf_vv {
            let interf = Interference::VertexVertex { v1: inf.v1, v2: inf.v2, merged_vertex: inf.merged_vertex };
            if Self::is_non_trivial(&interf, &ds, a_vc, a_ec, a_fc) && Self::is_allowed_by_level(&interf, level) {
                filtered.push(interf);
            }
        }
        for inf in &ds.interf_ve {
            let interf = Interference::VertexEdge { vertex: inf.vertex, edge: inf.edge, param: inf.param };
            if Self::is_non_trivial(&interf, &ds, a_vc, a_ec, a_fc) && Self::is_allowed_by_level(&interf, level) {
                filtered.push(interf);
            }
        }
        for inf in &ds.interf_vf {
            let interf = Interference::VertexFace { vertex: inf.vertex, face: inf.face };
            if Self::is_non_trivial(&interf, &ds, a_vc, a_ec, a_fc) && Self::is_allowed_by_level(&interf, level) {
                filtered.push(interf);
            }
        }
        for inf in &ds.interf_ee {
            let interf = Interference::EdgeEdge { e1: inf.e1, e2: inf.e2, point: inf.point, param1: inf.param1, param2: inf.param2, new_vertex: inf.new_vertex };
            if Self::is_non_trivial(&interf, &ds, a_vc, a_ec, a_fc) && Self::is_allowed_by_level(&interf, level) {
                filtered.push(interf);
            }
        }
        for inf in &ds.interf_ef {
            let interf = Interference::EdgeFace { edge: inf.edge, face: inf.face, point: inf.point, edge_param: inf.edge_param, new_vertex: inf.new_vertex };
            if Self::is_non_trivial(&interf, &ds, a_vc, a_ec, a_fc) && Self::is_allowed_by_level(&interf, level) {
                filtered.push(interf);
            }
        }
        for inf in &ds.interf_ff {
            let interf = Interference::FaceFace { f1: inf.f1, f2: inf.f2, curves: inf.curves.clone(), points: inf.points.clone() };
            if Self::is_non_trivial(&interf, &ds, a_vc, a_ec, a_fc) && Self::is_allowed_by_level(&interf, level) {
                filtered.push(interf);
            }
        }

        self.interferences = filtered;
        self.has_interferences = !self.interferences.is_empty();
    }

    /// Returns true if any self-interferences were found in the last Perform().
    ///
    /// OCCT ref: BOPAlgo_CheckerSI::HasInterferences (BOPAlgo_CheckerSI.cxx L220-240).
    pub fn has_interferences(&self) -> bool {
        self.has_interferences
    }

    /// Get the list of self-interferences found in the last Perform().
    ///
    /// OCCT ref: BOPAlgo_CheckerSI::GetInterferences (BOPAlgo_CheckerSI.cxx L250-260).
    pub fn get_interferences(&self) -> &[Interference] {
        &self.interferences
    }

    // ── Filter helpers ──────────────────────────────────────────────────────────

    /// Returns `true` if the interference represents a true self-interference
    /// (not a trivial A-vs-B coincidence from loading the same shape twice).
    ///
    /// Trivial interferences occur when loading the same shape twice: a vertex
    /// from copy A compared against the same vertex from copy B always coincides,
    /// but that is not a "self-interference" — it is just the same entity.
    fn is_non_trivial(
        interf: &Interference,
        ds: &DS,
        a_vc: usize,
        a_ec: usize,
        a_fc: usize,
    ) -> bool {
        match interf {
            Interference::VertexVertex { v1, v2, .. } => {
                Self::same_type_filter(*v1, *v2, a_vc)
            }
            Interference::VertexEdge { vertex, edge, .. } => {
                Self::ve_filter(*vertex, *edge, ds, a_vc)
            }
            Interference::EdgeEdge { e1, e2, .. } => {
                Self::same_type_filter(*e1, *e2, a_ec)
            }
            Interference::VertexFace { vertex, face } => {
                Self::vf_filter(*vertex, *face, ds, a_vc)
            }
            Interference::EdgeFace { edge, face, .. } => {
                Self::ef_filter(*edge, *face, ds, a_ec)
            }
            Interference::FaceFace { f1, f2, .. } => {
                Self::same_type_filter(*f1, *f2, a_fc)
            }
        }
    }

    /// Filter for same-type interferences (VV, EE, FF).
    ///
    /// Returns `true` if the two indices represent distinct original entities.
    /// When one index is from the A-copy range and the other from the B-copy
    /// range at a matching offset, they are the same original entity → trivial.
    fn same_type_filter(idx1: usize, idx2: usize, count: usize) -> bool {
        let (a_idx, b_idx) = if idx1 < idx2 {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };
        // Same entity: one is in range [0, count) and the other in [count, 2*count)
        // with matching offset.
        let is_same = a_idx < count && b_idx >= count && b_idx < 2 * count && a_idx + count == b_idx;
        !is_same
    }

    /// Filter for VertexEdge: returns `true` if the vertex is NOT an endpoint
    /// of the edge (i.e., it is a true self-interference, not natural topology).
    ///
    /// Since PaveFiller iterates both A→B and B→A, `vertex` may be from the A
    /// copy and `edge` from the B copy, or vice versa.
    fn ve_filter(vertex: usize, edge: usize, ds: &DS, a_vc: usize) -> bool {
        if edge >= ds.edges.len() {
            return true;
        }
        let e = &ds.edges[edge];
        // Determine which copy the vertex belongs to and check accordingly.
        // The edge's start/end vertices are from one copy; the vertex may be
        // from the *other* copy (with offset).
        if vertex < a_vc {
            // Vertex is from A copy. Edge may be from B copy (vertex offset = a_vc).
            let b_start = e.start_vertex.wrapping_sub(a_vc);
            let b_end = e.end_vertex.wrapping_sub(a_vc);
            vertex != b_start && vertex != b_end
        } else if vertex >= a_vc && vertex < 2 * a_vc {
            // Vertex is from B copy. Edge may be from A copy (vertex offset = a_vc added back).
            let a_start = e.start_vertex.wrapping_add(a_vc);
            let a_end = e.end_vertex.wrapping_add(a_vc);
            vertex != a_start && vertex != a_end
        } else {
            // Vertex is from intersection (index >= 2*a_vc). This is a true interference.
            true
        }
    }

    /// Filter for VertexFace: returns `true` if the vertex is NOT a boundary
    /// vertex of the face (i.e., it is truly on the face interior).
    fn vf_filter(vertex: usize, face: usize, ds: &DS, a_vc: usize) -> bool {
        if face >= ds.faces.len() {
            return true;
        }
        let f = &ds.faces[face];
        if vertex < a_vc {
            // Vertex is from A copy. Face boundary_verts are from B copy (offset = a_vc).
            for &bv in &f.boundary_verts {
                if bv >= a_vc && vertex + a_vc == bv {
                    return false;
                }
            }
        } else if vertex >= a_vc && vertex < 2 * a_vc {
            // Vertex is from B copy. Face boundary_verts may be from A copy.
            for &bv in &f.boundary_verts {
                if bv < a_vc && vertex == bv + a_vc {
                    return false;
                }
            }
        }
        true
    }

    /// Filter for EdgeFace: returns `true` if the edge is NOT a boundary edge
    /// of the face (i.e., it is a true edge-face self-interference).
    fn ef_filter(edge: usize, face: usize, ds: &DS, a_ec: usize) -> bool {
        if face >= ds.faces.len() {
            return true;
        }
        let f = &ds.faces[face];
        if edge < a_ec {
            // Edge is from A copy. Face boundary_edges are from B copy (offset = a_ec).
            for &be in &f.boundary_edges {
                if be >= a_ec && edge + a_ec == be {
                    return false;
                }
            }
        } else if edge >= a_ec && edge < 2 * a_ec {
            // Edge is from B copy. Face boundary_edges may be from A copy.
            for &be in &f.boundary_edges {
                if be < a_ec && edge == be + a_ec {
                    return false;
                }
            }
        }
        true
    }

    /// Check if an interference type is allowed by the current LevelOfCheck.
    ///
    /// OCCT BOPAlgo_CheckerSI LevelOfCheck:
    /// - Level 0: VE, EE, EF, FF (basic entity checks)
    /// - Level 1: +VV
    /// - Level 2: +VF
    /// - Level 3+: all
    fn is_allowed_by_level(interf: &Interference, level: u8) -> bool {
        match interf {
            Interference::VertexVertex { .. } => level >= 1,
            Interference::VertexFace { .. } => level >= 2,
            _ => true, // VE, EE, EF, FF are always included
        }
    }
}

impl Default for CheckerSI {
    fn default() -> Self {
        Self::new()
    }
}
