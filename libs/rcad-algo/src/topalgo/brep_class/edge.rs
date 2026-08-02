// OCCT BRepClass_Edge (BRepClass_Edge.hxx / .cxx)
// A boundary edge of a face, with the face it belongs to and the next edge
// sharing the last vertex (used by BRepClass_Intersector::CheckSkip).

use crate::topalgo::shape_source::ShapeSource;

/// OCCT BRepClass_Edge — edge + face references for 2D face classification.
pub struct ClassEdge {
    /// DS index of the edge.
    edge: usize,
    /// DS index of the face.
    face: usize,
    /// OCCT myMaxTolerance.
    max_tolerance: f64,
    /// OCCT myUseBndBox.
    use_bnd_box: bool,
    /// OCCT myNextEdge — the next edge in the wire sharing the last vertex.
    next_edge: Option<usize>,
}

impl ClassEdge {
    pub fn new(edge: usize, face: usize) -> Self {
        ClassEdge {
            edge,
            face,
            max_tolerance: f64::INFINITY,
            use_bnd_box: false,
            next_edge: None,
        }
    }

    pub fn edge(&self) -> usize {
        self.edge
    }

    pub fn face(&self) -> usize {
        self.face
    }

    pub fn max_tolerance(&self) -> f64 {
        self.max_tolerance
    }

    pub fn set_max_tolerance(&mut self, t: f64) {
        self.max_tolerance = t;
    }

    pub fn use_bnd_box(&self) -> bool {
        self.use_bnd_box
    }

    pub fn set_use_bnd_box(&mut self, b: bool) {
        self.use_bnd_box = b;
    }

    pub fn next_edge(&self) -> Option<usize> {
        self.next_edge
    }

    /// OCCT BRepClass_Edge::SetNextEdge (BRepClass_Edge.cxx L34-61) — find the
    /// next edge that shares the edge's last vertex. Only set when exactly two
    /// edges share that vertex.
    pub fn set_next_edge(&mut self, ds: &dyn ShapeSource) {
        if self.next_edge.is_some() || self.edge >= ds.nb_shapes() {
            return;
        }
        // The edge's sub-shapes are its two vertices in DS index order
        // [first, last] (BOPDS_DS::sub_shapes_of on an edge).
        let vsubs = ds.sub_shapes(self.edge);
        if vsubs.len() < 2 {
            return;
        }
        let a_vl = vsubs[vsubs.len() - 1];
        if a_vl == vsubs[0] {
            return; // closed/degenerate edge: same vertex at both ends
        }
        let edges_at_vl = ds.map_ve(a_vl).cloned().unwrap_or_default();
        if edges_at_vl.len() == 2 {
            for &other in &edges_at_vl {
                if other != self.edge {
                    self.next_edge = Some(other);
                    return;
                }
            }
        }
    }
}
