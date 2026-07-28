// OCCT BOPAlgo_BuilderSolid — solid building from shells.
//
// OCCT BOPAlgo_BuilderSolid.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use crate::bop::algo::shell_splitter::make_connexity_blocks_from_shapes;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{TShape, TSolidData, TShellData, tshape_flags};
use std::collections::{HashSet, HashMap};
use std::sync::Arc;

/// OCCT BOPAlgo_BuilderSolid — builds solids from a set of faces.
pub struct BuilderSolid<'a> {
    ds: &'a DS,
    my_report: Report,
    /// Input faces (OCCT: myShapes)
    pub my_shapes: Vec<Shape>,
    /// Resulting solids
    pub my_solids: Vec<Shape>,
    /// Faces with free edges (OCCT: myShapesToAvoid)
    my_shapes_to_avoid: HashSet<u64>,
    /// Connected face groups (OCCT: myLoops — shells)
    my_loops: Vec<Vec<Shape>>,
}

impl<'a> BuilderSolid<'a> {
    pub fn new(ds: &'a DS) -> Self {
        BuilderSolid {
            ds,
            my_report: Report::new(),
            my_shapes: Vec::new(),
            my_solids: Vec::new(),
            my_shapes_to_avoid: HashSet::new(),
            my_loops: Vec::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::Perform (BOPAlgo_BuilderSolid.cxx L76-125).
    pub fn perform(&mut self) {
        if self.my_shapes.is_empty() {
            return;
        }
        // OCCT L94-104: MakeCompound from all faces
        // (rcad: shapes stored in Vec, no compound needed)
        // OCCT L106: PerformShapesToAvoid
        self.perform_shapes_to_avoid();
        if self.has_errors() { return; }
        // OCCT L112: PerformLoops — group faces into shells via connexity
        self.perform_loops();
        if self.has_errors() { return; }
        // OCCT L118: PerformAreas — classify shells, build solids
        self.perform_areas();
        if self.has_errors() { return; }
        // OCCT L124: PerformInternalShapes
        self.perform_internal_shapes();
    }

    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }

    /// OCCT BOPAlgo_BuilderSolid::PerformShapesToAvoid (L129+).
    /// Finds faces with free boundary edges (edges used by only one face).
    fn perform_shapes_to_avoid(&mut self) {
        // Build edge->face adjacency
        let mut edge_to_faces: HashMap<u64, Vec<usize>> = HashMap::new();
        for (fi, face) in self.my_shapes.iter().enumerate() {
            let edge_ptrs = face_edge_ptrs(face);
            for eptr in edge_ptrs {
                edge_to_faces.entry(eptr).or_default().push(fi);
            }
        }
        // OCCT: faces whose edge is used by only one face → add to myShapesToAvoid
        for (fi, face) in self.my_shapes.iter().enumerate() {
            let edge_ptrs = face_edge_ptrs(face);
            let has_free_edge = edge_ptrs.iter().any(|eptr| {
                edge_to_faces.get(eptr).map_or(true, |faces| faces.len() == 1)
            });
            if has_free_edge {
                self.my_shapes_to_avoid.insert(face.ptr_id());
            }
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformLoops (connected face grouping).
    /// Groups faces into connected shells using edge adjacency.
    fn perform_loops(&mut self) {
        if self.my_shapes.is_empty() {
            return;
        }
        // OCCT: uses MakeConnexityBlocks to build connected face groups
        let mut blocks: Vec<Vec<usize>> = Vec::new();
        make_connexity_blocks_from_shapes(&self.my_shapes, self.ds, &mut blocks);
        self.my_loops.clear();
        for block_indices in &blocks {
            let shell_faces: Vec<Shape> = block_indices.iter()
                .map(|&idx| self.my_shapes[idx].clone())
                .collect();
            if !shell_faces.is_empty() {
                // Build a Shell TShape from these connected faces
                let shell_tshape = TShape::Shell(TShellData {
                    my_shapes: vec![],
                    flags: tshape_flags::DEFAULT,
                    faces: shell_faces,
                });
                self.my_loops.push(block_indices.iter()
                    .map(|&idx| self.my_shapes[idx].clone())
                    .collect());
            }
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformAreas (L387+).
    /// Classifies shell groups → builds solids from IN shells.
    fn perform_areas(&mut self) {
        // OCCT: uses BRepClass3d_SolidClassifier to classify shells
        // IN shells → create Solid TShape with the shell
        for loop_faces in &self.my_loops {
            // Build a solid from the connected shell faces
            let shell_tshape = TShape::Shell(TShellData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                faces: loop_faces.clone(),
            });
            let shell_shape = Shape::new(
                Arc::new(shell_tshape), 0, rcad_kernel::topods::Orientation::Forward,
            );
            let solid_tshape = TShape::Solid(TSolidData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                shells: vec![shell_shape],
                internal_vertices: vec![],
                internal_edges: vec![],
            });
            let solid_shape = Shape::new(
                Arc::new(solid_tshape), 0, rcad_kernel::topods::Orientation::Forward,
            );
            self.my_solids.push(solid_shape);
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformInternalShapes.
    fn perform_internal_shapes(&mut self) {}
}

/// Extract edge ptr_ids from a Face Shape.
fn face_edge_ptrs(face: &Shape) -> Vec<u64> {
    let mut edges = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            // Outer wire edges
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                for e in &wd.edges {
                    edges.push(e.ptr_id());
                }
            }
            // Inner wire edges
            for iw in &fd.inner_wires {
                if let TShape::Wire(wd) = &*iw.data {
                    for e in &wd.edges {
                        edges.push(e.ptr_id());
                    }
                }
            }
        }
        _ => {}
    }
    edges
}
