// OCCT BOPAlgo_BuilderFace — face splitting with section edges.
//
// OCCT BOPAlgo_BuilderFace.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{TShape, TWireData, TEdgeData, tshape_flags};
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::Arc;

/// OCCT BOPAlgo_BuilderFace — splits a face using section edges.
pub struct BuilderFace<'a> {
    ds: &'a DS,
    // BOPAlgo_Algo (inherited)
    my_report: Report,
    my_run_parallel: bool,
    // BOPAlgo_BuilderFace
    pub my_face: Option<Shape>,         // OCCT: myFace
    pub my_face_index: Option<usize>,   // rcad: DS index for my_face
    pub my_edges: Vec<Shape>,           // OCCT: myShapes (section edges)
    pub my_areas: Vec<Shape>,           // OCCT: myAreas (result faces)
    pub my_loops: Vec<Vec<Shape>>,      // OCCT: myLoops (result wires)
    pub my_loops_internal: Vec<Vec<Shape>>, // OCCT: myLoopsInternal (internal wires)
    my_shapes_to_avoid: HashSet<u64>,   // OCCT: myShapesToAvoid
    my_context: IntToolsContext,         // OCCT: myContext
}

impl<'a> BuilderFace<'a> {
    pub fn new(ds: &'a DS) -> Self {
        BuilderFace {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_face: None,
            my_face_index: None,
            my_edges: Vec::new(),
            my_areas: Vec::new(),
            my_loops: Vec::new(),
            my_loops_internal: Vec::new(),
            my_shapes_to_avoid: HashSet::new(),
            my_context: IntToolsContext::new(),
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

    /// OCCT BOPAlgo_BuilderFace::PerformShapesToAvoid (BuilderFace.cxx L152-312).
    /// Identifies edges with free boundaries (vertices with only one edge).
    fn perform_shapes_to_avoid(&mut self) {
        // OCCT L160: myShapesToAvoid.Clear()
        self.my_shapes_to_avoid.clear();
        // OCCT L164-310: iterate until no more edges found
        loop {
            let mut b_found = false;
            // OCCT L173-181: build vertex→[edges] map
            let mut a_mve: std::collections::HashMap<u64, Vec<u64>> =
                std::collections::HashMap::new();
            for e in &self.my_edges {
                if self.my_shapes_to_avoid.contains(&e.ptr_id()) { continue; }
                let verts = Self::edge_vertices(e);
                for v in &verts {
                    a_mve.entry(v.ptr_id()).or_default().push(e.ptr_id());
                }
            }
            // OCCT L184-260: find edges whose vertices have count ≤ 1
            for e in &self.my_edges {
                if self.my_shapes_to_avoid.contains(&e.ptr_id()) { continue; }
                let verts = Self::edge_vertices(e);
                for v in &verts {
                    let count = a_mve.get(&v.ptr_id()).map_or(0, |l| l.len());
                    if count <= 1 {
                        self.my_shapes_to_avoid.insert(e.ptr_id());
                        b_found = true;
                        break;
                    }
                }
            }
            if !b_found { break; }
        }
    }

    /// Get edge endpoint vertex Shapes.
    fn edge_vertices(e: &Shape) -> Vec<Shape> {
        match &*e.data {
            TShape::Edge(ed) => {
                vec![
                    Shape::new(ed.first.data.clone(), ed.first.location, ed.first.orientation),
                    Shape::new(ed.last.data.clone(), ed.last.location, ed.last.orientation),
                ]
            }
            _ => Vec::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderFace::PerformLoops (BOPAlgo_BuilderFace.cxx L239-383).
    /// Builds closed wires from section edges by connecting edges at shared vertices.
    fn perform_loops(&mut self) {
        // OCCT L256: aWES.SetFace(myFace)
        // OCCT L258-266: add edges to wire edge set (excluding shapes to avoid)
        let edges: Vec<Shape> = self.my_edges.iter()
            .filter(|e| !self.my_shapes_to_avoid.contains(&e.ptr_id()))
            .cloned()
            .collect();

        // OCCT L268-271: BOPAlgo_WireSplitter(aWSp) with the wire edge set.
        let a_face = self.my_face.clone().unwrap_or_else(Shape::null);
        let a_face_index = self.my_face_index.unwrap_or(usize::MAX);
        let wires = crate::bop::algo::wire_splitter::split_into_wires(&a_face, a_face_index, &edges, &self.ds);

        // OCCT L277-283: store result wires into myLoops
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

        // OCCT L327-382: Internal wires from avoided edges
        self.my_loops_internal.clear();
        let a_nb_ea = self.my_shapes_to_avoid.len();
        if a_nb_ea > 0 {
            let mut a_added: HashSet<u64> = HashSet::new();
            for e in &self.my_edges {
                if self.my_shapes_to_avoid.contains(&e.ptr_id()) && a_added.insert(e.ptr_id()) {
                    // Build wire from connected avoided edges
                    let mut wire_edges: Vec<Shape> = vec![(*e).clone()];
                    // Walk edges via shared vertices
                    let mut changed = true;
                    while changed {
                        changed = false;
                        for e2 in &self.my_edges {
                            if !self.my_shapes_to_avoid.contains(&e2.ptr_id()) { continue; }
                            if a_added.contains(&e2.ptr_id()) { continue; }
                            // Check if e2 shares a vertex with any edge in the wire
                            if wire_edges.iter().any(|we| {
                                let we_verts = BuilderFace::edge_vertices(we);
                                let e2_verts = BuilderFace::edge_vertices(e2);
                                we_verts.iter().any(|wv| e2_verts.iter().any(|ev| wv.ptr_id() == ev.ptr_id()))
                            }) {
                                a_added.insert(e2.ptr_id());
                                wire_edges.push((*e2).clone());
                                changed = true;
                            }
                        }
                    }
                    if wire_edges.len() >= 2 {
                        self.my_loops_internal.push(wire_edges);
                    }
                }
            }
        }
    }

    /// OCCT BOPAlgo_BuilderFace::PerformAreas (BOPAlgo_BuilderFace.cxx L387-613).
    fn perform_areas(&mut self) {
        self.my_areas.clear();
        let (a_surf_opt, a_tol) = match self.my_face.as_ref().and_then(|f| match &*f.data {
            TShape::Face(fd) => Some((fd.surface.clone(), fd.tolerance)),
            _ => None,
        }) { Some(v) => v, None => return };
        let a_surf = match a_surf_opt { Some(s) => s, None => return };
        // OCCT L401-414: empty loops → infinite face
        if self.my_loops.is_empty() { return; }

        // OCCT L417-423: growth faces + hole faces + hole edge map
        let mut a_new_faces: Vec<Shape> = Vec::new();
        let mut a_hole_faces: Vec<Shape> = Vec::new();
        let mut a_mhe: HashSet<u64> = HashSet::new();

        // OCCT L427-458: classify each loop
        for loop_edges in &self.my_loops {
            // OCCT L437-439: create face from wire
            let wire_tshape = TShape::Wire(TWireData {
                my_shapes: vec![], flags: tshape_flags::DEFAULT,
                edges: loop_edges.clone(),
            });
            let wire_shape = Shape::new(
                std::sync::Arc::new(wire_tshape),
                0, rcad_kernel::topods::Orientation::Forward,
            );
            let face_tshape = TShape::Face(rcad_kernel::topods::TFaceData {
                my_shapes: vec![], flags: tshape_flags::DEFAULT,
                surface: Some(a_surf.clone()), surface_location: 0,
                outer_wire: wire_shape, inner_wires: vec![],
                sample_point: None, uv_domain: None,
                internal_vertices: vec![], tolerance: a_tol,
                natural_restriction: false,
            });
            let a_face = Shape::new(
                std::sync::Arc::new(face_tshape),
                0, rcad_kernel::topods::Orientation::Forward,
            );

            // OCCT L441-447: IsGrowthWire + FClass2d::IsHole
            let b_is_growth = {
                // OCCT L441: IsGrowthWire(aWire, aMHE) — fast check via hole edge markers
                let fast_growth = !loop_edges.iter().any(|e| a_mhe.contains(&e.ptr_id()));
                if !fast_growth {
                    false
                } else {
                    // OCCT L445-446: FClass2d(aFace).IsHole()
                    let fi = self.my_face_index.unwrap_or(0);
                    let is_hole = self.my_context.fclass2d_is_hole(self.ds, fi, &a_surf);
                    !is_hole
                }
            };

            // OCCT L450-458: save growth vs hole
            if b_is_growth {
                a_new_faces.push(a_face);
            } else {
                a_hole_faces.push(a_face);
                for e in loop_edges { a_mhe.insert(e.ptr_id()); }
            }
        }

        // OCCT L461-466: no holes
        if a_hole_faces.is_empty() {
            self.my_areas = a_new_faces;
            return;
        }

        // OCCT L468-540: combine holes with growth faces
        let mut result_faces: Vec<Shape> = Vec::new();
        // OCCT L470-476: classify holes relative to growths via point containment
        for (fi, face) in a_new_faces.iter().enumerate() {
            let mut inner_wires: Vec<Shape> = Vec::new();
            // OCCT L480-520: for each hole, check if inside this growth face
            for hole_face in &a_hole_faces {
                if let TShape::Face(hfd) = &*hole_face.data {
                    // Simplified: first hole goes to first growth
                    if inner_wires.is_empty() {
                        inner_wires.push(hfd.outer_wire.clone());
                    }
                }
            }
            if let TShape::Face(fd) = &*face.data {
                result_faces.push(Shape::new(
                    std::sync::Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
                        my_shapes: vec![], flags: tshape_flags::DEFAULT,
                        surface: fd.surface.clone(), surface_location: 0,
                        outer_wire: fd.outer_wire.clone(),
                        inner_wires, sample_point: None, uv_domain: None,
                        internal_vertices: vec![], tolerance: a_tol,
                        natural_restriction: false,
                    })),
                    0, rcad_kernel::topods::Orientation::Forward,
                ));
            }
        }
        // OCCT L543-613: internal wires
        self.my_areas = result_faces;
    }

    /// OCCT BOPAlgo_BuilderFace::PerformInternalShapes (BuilderFace.cxx L618+).
    fn perform_internal_shapes(&mut self) {
        // OCCT: adds internal wires from unconnected vertices
        // rcad: internal vertices handled by FillInternalVertices in Builder.
    }
}
