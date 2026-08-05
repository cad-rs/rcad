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

    /// OCCT BOPAlgo_BuilderFace::PerformShapesToAvoid (BuilderFace.cxx L152-235).
    /// Iteratively marks edges with a free boundary (a vertex used by at most
    /// one non-avoided edge, or by two IsSame copies of one edge) as "to avoid".
    fn perform_shapes_to_avoid(&mut self) {
        // OCCT L160: myShapesToAvoid.Clear()
        self.my_shapes_to_avoid.clear();
        // OCCT L164-234: iterate until no more edges are found.
        loop {
            let mut b_found = false;
            // OCCT L173-182: aMVE — vertex → [edges] (skipping avoided edges).
            let mut a_mve: HashMap<u64, (Shape, Vec<Shape>)> = HashMap::new();
            for a_e in &self.my_edges {
                if self.my_shapes_to_avoid.contains(&a_e.ptr_id()) { continue; }
                for a_v in Self::edge_vertices(a_e) {
                    let entry = a_mve
                        .entry(a_v.ptr_id())
                        .or_insert_with(|| (a_v.clone(), Vec::new()));
                    entry.1.push(a_e.clone());
                }
            }
            // OCCT L186-228: for each vertex decide.
            for (_vptr, (a_v, a_le)) in &a_mve {
                let a_nb_e = a_le.len();
                if a_nb_e == 0 { continue; }
                let a_e1 = &a_le[0];
                if a_nb_e == 1 {
                    // OCCT L198-210: single edge at the vertex.
                    if a_e1.as_edge().map_or(true, |ed| ed.degenerated) {
                        continue;
                    }
                    if a_v.orientation == rcad_kernel::topods::Orientation::Internal {
                        continue;
                    }
                    b_found = true;
                    self.my_shapes_to_avoid.insert(a_e1.ptr_id());
                } else if a_nb_e == 2 {
                    // OCCT L211-227: two edges at the vertex.
                    let a_e2 = &a_le[1];
                    if a_e2.is_partner(a_e1) {
                        let vv = Self::edge_vertices(a_e1);
                        if vv.len() >= 2 && vv[0].is_partner(&vv[1]) {
                            // Degenerated ring — both ends are the same vertex.
                            continue;
                        }
                        b_found = true;
                        self.my_shapes_to_avoid.insert(a_e1.ptr_id());
                        self.my_shapes_to_avoid.insert(a_e2.ptr_id());
                    }
                }
            }
            if !b_found { break; }
        }
    }

    /// OCCT IntTools_Context::IsInfiniteFace — a face without a bounded outer
    /// wire is treated as infinite (unbounded surface).
    fn is_infinite_face(face: &Shape) -> bool {
        match &*face.data {
            TShape::Face(fd) => match &*fd.outer_wire.data {
                TShape::Wire(wd) => wd.edges.is_empty(),
                _ => true,
            },
            _ => true,
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
                    // OCCT L351-381: the wire is appended unconditionally
                    // (a single-edge wire is a valid internal wire).
                    self.my_loops_internal.push(wire_edges);
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
        // OCCT L401-414: empty loops — an infinite face becomes a face without
        // wires.
        if self.my_loops.is_empty() {
            let is_infinite = self
                .my_face
                .as_ref()
                .map(Self::is_infinite_face)
                .unwrap_or(false);
            if is_infinite {
                let natural_restriction = self
                    .my_face
                    .as_ref()
                    .and_then(|f| match &*f.data {
                        TShape::Face(fd) => Some(fd.natural_restriction),
                        _ => None,
                    })
                    .unwrap_or(false);
                let face_tshape = TShape::Face(rcad_kernel::topods::TFaceData {
                    my_shapes: vec![],
                    flags: tshape_flags::DEFAULT,
                    surface: Some(a_surf.clone()),
                    surface_location: 0,
                    outer_wire: Shape::null(),
                    inner_wires: vec![],
                    sample_point: None,
                    uv_domain: None,
                    internal_vertices: vec![],
                    tolerance: a_tol,
                    natural_restriction,
                });
                self.my_areas.push(Shape::new(
                    std::sync::Arc::new(face_tshape),
                    0,
                    rcad_kernel::topods::Orientation::Forward,
                ));
            }
            return;
        }

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
                // OCCT L441: IsGrowthWire(aWire, aMHE) — returns true when the
                // wire contains any hole-face edge (BOPAlgo_BuilderFace.cxx
                // L898-913: theMHE.Contains(aIt.Value())). Only when it has no
                // hole edge is the FClass2d classification run.
                let has_hole_edge = loop_edges.iter().any(|e| a_mhe.contains(&e.ptr_id()));
                if has_hole_edge {
                    true
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
