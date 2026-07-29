// OCCT BRepClass_FaceExplorer (BRepClass_FaceExplorer.hxx / .cxx)
// Provides exploration of a face's 2D edges for classification.

use glam::DVec2;

/// OCCT BRepClass_Edge — 2D edge representation.
pub struct ClassEdge {
    pub start: DVec2,
    pub end: DVec2,
}

/// OCCT BRepClass_FaceExplorer — explores a face's wires and edges in 2D.
///
/// rcad: simplified — provides edge data from face pcurves.
pub struct FaceExplorer {
    // Index into wires/edges
    current_wire: usize,
    current_edge: usize,
    // Edge data
    edges: Vec<Vec<ClassEdge>>,
    // Tolerance
    max_tolerance: f64,
}

impl FaceExplorer {
    /// OCCT: Constructor(F) — initialize from face.
    pub fn new() -> Self {
        FaceExplorer {
            current_wire: 0, current_edge: 0,
            edges: Vec::new(), max_tolerance: 1e-7,
        }
    }

    /// Initialize with face data.
    pub fn init(&mut self, edges_per_wire: Vec<Vec<ClassEdge>>) {
        self.edges = edges_per_wire;
        self.current_wire = 0;
        self.current_edge = 0;
    }

    /// OCCT: CheckPoint(P) — adjust point if too far from bounding box.
    pub fn check_point(&self, _point: &mut DVec2) -> bool { true }

    /// OCCT: Reject(P) — quick bounding box rejection.
    pub fn reject(&self, _p: DVec2) -> bool { false }

    /// OCCT: Segment(P, L, Par) — build a 2D segment for intersection.
    pub fn segment(&self, _p: DVec2) -> Option<(glam::DVec2, glam::DVec2, f64)> { None }

    /// OCCT: OtherSegment(P, L, Par) — alternative segment.
    pub fn other_segment(&self, _p: DVec2) -> Option<(glam::DVec2, glam::DVec2, f64)> { None }

    /// Wire iteration.
    pub fn init_wires(&mut self) { self.current_wire = 0; }
    pub fn more_wires(&self) -> bool { self.current_wire < self.edges.len() }
    pub fn next_wire(&mut self) { self.current_wire += 1; }
    pub fn reject_wire(&self, _l: glam::DVec2, _par: f64) -> bool { false }

    /// Edge iteration within current wire.
    pub fn init_edges(&mut self) { self.current_edge = 0; }
    pub fn more_edges(&self) -> bool {
        if self.current_wire < self.edges.len() {
            self.current_edge < self.edges[self.current_wire].len()
        } else { false }
    }
    pub fn next_edge(&mut self) { self.current_edge += 1; }
    pub fn reject_edge(&self, _l: glam::DVec2, _par: f64) -> bool { false }

    /// Current edge.
    pub fn current_edge_data(&self) -> Option<(&ClassEdge, u8)> {
        if self.current_wire < self.edges.len()
            && self.current_edge < self.edges[self.current_wire].len() {
            Some((&self.edges[self.current_wire][self.current_edge], 0))
        } else { None }
    }

    pub fn max_tolerance(&self) -> f64 { self.max_tolerance }
    pub fn set_max_tolerance(&mut self, tol: f64) { self.max_tolerance = tol; }
}
