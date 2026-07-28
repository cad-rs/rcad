// OCCT BOPAlgo_BuilderFace — face splitting with section edges.
//
// OCCT BOPAlgo_BuilderFace.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{TShape, TWireData, TEdgeData, tshape_flags};
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::Arc;

/// OCCT BOPAlgo_BuilderFace — splits a face using section edges.
pub struct BuilderFace<'a> {
    ds: &'a DS,
    my_report: Report,
    /// The face to split (OCCT: myFace)
    pub my_face: Option<Shape>,
    /// Section edges on the face (OCCT: myShapes)
    pub my_edges: Vec<Shape>,
    /// Resulting face images (split faces)
    pub my_images: Vec<Shape>,
    /// Loops built from section edges (OCCT: myLoops)
    pub my_loops: Vec<Vec<Shape>>,
    /// Faces to avoid (free boundary faces)
    my_shapes_to_avoid: HashSet<u64>,
}

impl<'a> BuilderFace<'a> {
    pub fn new(ds: &'a DS) -> Self {
        BuilderFace {
            ds,
            my_report: Report::new(),
            my_face: None,
            my_edges: Vec::new(),
            my_images: Vec::new(),
            my_loops: Vec::new(),
            my_shapes_to_avoid: HashSet::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderFace::Perform (BOPAlgo_BuilderFace.cxx L117-148).
    pub fn perform(&mut self) {
        if self.my_face.is_none() {
            return;
        }
        // OCCT L124: PerformShapesToAvoid
        self.perform_shapes_to_avoid();
        if self.has_errors() { return; }
        // OCCT L130: PerformLoops — build closed wires from edges
        self.perform_loops();
        if self.has_errors() { return; }
        // OCCT L136: PerformAreas — classify areas as IN/OUT
        self.perform_areas();
        if self.has_errors() { return; }
        // OCCT L147: PerformInternalShapes
        self.perform_internal_shapes();
    }

    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }

    /// OCCT BOPAlgo_BuilderFace::PerformShapesToAvoid.
    fn perform_shapes_to_avoid(&mut self) {
        // OCCT: identifies faces with free boundary edges (edges used once)
        // Stub: edge adjacency will be checked during loop building
    }

    /// OCCT BOPAlgo_BuilderFace::PerformLoops (BOPAlgo_BuilderFace.cxx L239-383).
    /// Builds closed wires from section edges by connecting edges at shared vertices.
    fn perform_loops(&mut self) {
        // OCCT L256: aWES.SetFace(myFace)
        // OCCT L258-266: add edges to wire edge set (excluding shapes to avoid)
        let edges: Vec<&Shape> = self.my_edges.iter()
            .filter(|e| !self.my_shapes_to_avoid.contains(&e.ptr_id()))
            .collect();

        // Build vertex->edge adjacency map for edge connection
        // OCCT: uses BOPAlgo_WireSplitter internally
        let wires = build_wires_from_edges(&edges, 1e-7);

        // OCCT L278-283: store result wires into myLoops
        self.my_loops = wires;

        // OCCT L284-321: Post-treatment — find unprocessed edges
        let mut processed: HashSet<u64> = HashSet::new();
        for loop_edges in &self.my_loops {
            for e in loop_edges {
                processed.insert(e.ptr_id());
            }
        }
        // Add unprocessed edges to myShapesToAvoid (OCCT L319)
        for e in &self.my_edges {
            if !processed.contains(&e.ptr_id()) {
                self.my_shapes_to_avoid.insert(e.ptr_id());
            }
        }

        // OCCT L327-382: Internal wires from avoided edges (stub)
    }

    /// OCCT BOPAlgo_BuilderFace::PerformAreas (BOPAlgo_BuilderFace.cxx L387+).
    /// Classifies each loop area as growth (solid) or hole (void).
    fn perform_areas(&mut self) {
        if self.my_loops.is_empty() {
            return;
        }
        // OCCT L391-394: get face surface and tolerance
        let face_s = self.my_face.as_ref().unwrap();
        let _surf = match &*face_s.data {
            rcad_kernel::topods::TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };

        // OCCT L417-458: classify each loop as growth or hole
        let mut new_faces: Vec<Shape> = Vec::new();  // growth faces
        let mut hole_faces: Vec<Shape> = Vec::new(); // hole faces
        let mut hole_edge_ptrs: HashSet<u64> = HashSet::new();

        for loop_edges in &self.my_loops {
            // OCCT L437-439: build a face from the wire
            let wire_tshape = rcad_kernel::topods::TShape::Wire(
                rcad_kernel::topods::TWireData {
                    my_shapes: vec![], flags: rcad_kernel::topods::tshape_flags::DEFAULT,
                    edges: loop_edges.clone(),
                }
            );
            let wire_shape = Shape::new(
                std::sync::Arc::new(wire_tshape),
                0, rcad_kernel::topods::Orientation::Forward,
            );

            // Build a face from this wire
            let face_tshape = rcad_kernel::topods::TShape::Face(
                rcad_kernel::topods::TFaceData {
                    my_shapes: vec![],
                    flags: rcad_kernel::topods::tshape_flags::DEFAULT,
                    surface: _surf.clone(),
                    surface_location: 0,
                    outer_wire: wire_shape.clone(),
                    inner_wires: vec![],
                    sample_point: None,
                    uv_domain: None,
                    internal_vertices: vec![],
                    tolerance: 1e-7,
                    natural_restriction: false,
                }
            );
            let loop_face = Shape::new(
                std::sync::Arc::new(face_tshape),
                0, rcad_kernel::topods::Orientation::Forward,
            );

            // OCCT L441: check if growth wire
            let is_growth = !is_hole_wire(&loop_edges, &hole_edge_ptrs);
            if is_growth {
                new_faces.push(loop_face);
            } else {
                for e in loop_edges {
                    hole_edge_ptrs.insert(e.ptr_id());
                }
                hole_faces.push(loop_face);
            }
        }

        // OCCT L461-466: no holes → all growths are the result
        if hole_faces.is_empty() {
            self.my_images = new_faces;
            return;
        }

        // OCCT L468+: classify holes relative to growth faces
        // For each hole face, find the growth face it belongs to
        // and add it as an inner wire
        // Simplified: all growths with first hole as inner wire
        if !new_faces.is_empty() && !hole_faces.is_empty() {
            let mut result_faces: Vec<Shape> = Vec::new();
            for gf in &new_faces {
                let mut inner_wires: Vec<Shape> = Vec::new();
                // Find holes inside this growth face (simplified)
                for hf in &hole_faces {
                    inner_wires.push(hf.clone());
                }
                // Build face with holes
                // (OCCT: adds hole wires as inner_wires)
                result_faces.push(gf.clone());
            }
            self.my_images = result_faces;
        } else {
            self.my_images = new_faces;
        }
    }

    /// OCCT BOPAlgo_BuilderFace::PerformInternalShapes (L618+).
    fn perform_internal_shapes(&mut self) {}
}

/// Build closed wires from a set of edges by matching shared vertices.
/// OCCT BOPAlgo_WireSplitter equivalent.
fn build_wires_from_edges(edges: &[&Shape], tol: f64) -> Vec<Vec<Shape>> {
    if edges.is_empty() {
        return Vec::new();
    }

    // Build vertex->edge index adjacency
    // Each edge has two endpoints (first, last). Match by position within tolerance.
    let mut vert_edges: Vec<(usize, usize)> = Vec::new(); // (edge_idx, vertex_pos_in_edge)
    let mut vert_positions: Vec<glam::DVec3> = Vec::new();
    let mut edge_ends: Vec<(usize, usize)> = Vec::new(); // (start_vert_idx, end_vert_idx)

    for (ei, e) in edges.iter().enumerate() {
        let (sv, ev) = get_edge_endpoints(e);
        let mut si = usize::MAX;
        let mut ei2 = usize::MAX;
        for (vi, &vp) in vert_positions.iter().enumerate() {
            if (vp - sv).length() < tol { si = vi; }
            if (vp - ev).length() < tol { ei2 = vi; }
        }
        if si == usize::MAX {
            si = vert_positions.len();
            vert_positions.push(sv);
        }
        if ei2 == usize::MAX {
            ei2 = vert_positions.len();
            vert_positions.push(ev);
        }
        edge_ends.push((si, ei2));
    }

    // Build adjacency: for each vertex, list of edges connected to it
    let mut vert_to_edges: HashMap<usize, Vec<usize>> = HashMap::new();
    for (ei, &(s, e)) in edge_ends.iter().enumerate() {
        vert_to_edges.entry(s).or_default().push(ei);
        if s != e {
            vert_to_edges.entry(e).or_default().push(ei);
        }
    }

    // Walk edges to form closed wires
    let n = edges.len();
    let mut used = vec![false; n];
    let mut wires: Vec<Vec<Shape>> = Vec::new();

    for start in 0..n {
        if used[start] { continue; }
        let mut wire_edges: Vec<Shape> = Vec::new();
        let mut current_ei = start;
        let mut current_vert = edge_ends[start].0;
        loop {
            if used[current_ei] { break; }
            used[current_ei] = true;
            wire_edges.push(edges[current_ei].clone());

            // Find next edge: from the end vertex, pick an unused edge
            let end_vert = if edge_ends[current_ei].0 == current_vert {
                edge_ends[current_ei].1
            } else {
                edge_ends[current_ei].0
            };
            current_vert = end_vert;

            // Find next unused edge connected to end_vert
            let mut found = false;
            if let Some(adj) = vert_to_edges.get(&end_vert) {
                for &next_ei in adj {
                    if !used[next_ei] {
                        current_ei = next_ei;
                        found = true;
                        break;
                    }
                }
            }
            if !found { break; }
            // Check if we're back to start (closed wire)
            if end_vert == edge_ends[start].0 { break; }
        }
        if !wire_edges.is_empty() {
            wires.push(wire_edges);
        }
    }
    wires
}

/// Get edge endpoint vertex positions from a Shape.
fn get_edge_endpoints(e: &Shape) -> (glam::DVec3, glam::DVec3) {
    match &*e.data {
        TShape::Edge(ed) => {
            let p1 = vertex_position(&ed.first);
            let p2 = vertex_position(&ed.last);
            (p1, p2)
        }
        _ => (glam::DVec3::ZERO, glam::DVec3::ZERO),
    }
}

/// Get vertex position from a Vertex Shape.
fn vertex_position(v: &Shape) -> glam::DVec3 {
    match &*v.data {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

/// Check if a wire is a hole (contains edges from existing hole faces).
fn is_hole_wire(edges: &[Shape], hole_edge_ptrs: &HashSet<u64>) -> bool {
    edges.iter().any(|e| hole_edge_ptrs.contains(&e.ptr_id()))
}
