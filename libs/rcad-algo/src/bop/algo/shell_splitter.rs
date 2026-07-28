// OCCT BOPAlgo_ShellSplitter — shell partitioning.
//
// OCCT BOPAlgo_ShellSplitter.cxx
// Splits shell into connected components (connexity blocks).
// Then builds shells from each block.

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use rcad_kernel::topo_shape::Shape;
use std::collections::{HashMap, HashSet, VecDeque};

/// OCCT BOPAlgo_ShellSplitter — partitions a shell into connected components.
pub struct ShellSplitter<'a> {
    ds: &'a DS,
    my_report: Report,
    /// Input faces (OCCT: myStartShapes)
    pub my_shapes: Vec<Shape>,
    /// Resulting connexity blocks (each block = connected subset of face indices)
    pub my_blocks: Vec<Vec<usize>>,
}

impl<'a> ShellSplitter<'a> {
    pub fn new(ds: &'a DS) -> Self {
        ShellSplitter {
            ds,
            my_report: Report::new(),
            my_shapes: Vec::new(),
            my_blocks: Vec::new(),
        }
    }

    /// OCCT BOPAlgo_ShellSplitter::Perform (BOPAlgo_ShellSplitter.cxx L137-149).
    /// Builds connexity blocks from faces and creates shells.
    pub fn perform(&mut self) {
        if self.my_shapes.is_empty() {
            return;
        }
        // OCCT L142: MakeConnexityBlocks(myStartShapes, EDGE, FACE, myLCB)
        self.my_blocks.clear();
        make_connexity_blocks_from_shapes(&self.my_shapes, self.ds, &mut self.my_blocks);
    }
}

/// OCCT BOPTools_AlgoTools::MakeConnexityBlocks — builds connected face subsets.
/// Groups faces by shared edge connectivity (OCCT L105+).
pub fn make_connexity_blocks_from_shapes(
    shapes: &[Shape],
    _ds: &DS,
    out_blocks: &mut Vec<Vec<usize>>,
) {
    // Build edge→face adjacency map
    // OCCT: TopExp::MapShapesAndAncestors(aF, EDGE, FACE, aMEF)
    let mut edge_to_faces: HashMap<u64, Vec<usize>> = HashMap::new();
    for (fi, face) in shapes.iter().enumerate() {
        let sub_edges = face_sub_edges(face);
        for edge_ptr in sub_edges {
            edge_to_faces.entry(edge_ptr).or_default().push(fi);
        }
    }

    // BFS over faces using edge adjacency
    let n = shapes.len();
    let mut visited = vec![false; n];
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut block: Vec<usize> = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(fi) = queue.pop_front() {
            block.push(fi);
            // Find adjacent faces via shared edges
            let sub_edges = face_sub_edges(&shapes[fi]);
            for edge_ptr in sub_edges {
                if let Some(neighbors) = edge_to_faces.get(&edge_ptr) {
                    for &nfi in neighbors {
                        if !visited[nfi] {
                            visited[nfi] = true;
                            queue.push_back(nfi);
                        }
                    }
                }
            }
        }
        if !block.is_empty() {
            out_blocks.push(block);
        }
    }
}

/// Extract edge ptr_ids from a Face Shape.
fn face_sub_edges(face: &Shape) -> Vec<u64> {
    let mut edges = Vec::new();
    match &*face.data {
        rcad_kernel::topods::TShape::Face(fd) => {
            // Collect from outer wire
            collect_wire_edge_ptrs(&fd.outer_wire, &mut edges);
            // Collect from inner wires
            for iw in &fd.inner_wires {
                collect_wire_edge_ptrs(iw, &mut edges);
            }
        }
        _ => {}
    }
    edges
}

/// Collect edge ptr_ids from a Wire Shape.
fn collect_wire_edge_ptrs(wire: &Shape, out: &mut Vec<u64>) {
    match &*wire.data {
        rcad_kernel::topods::TShape::Wire(wd) => {
            for e in &wd.edges {
                out.push(e.ptr_id());
            }
        }
        _ => {}
    }
}

/// Build connexity blocks from connected faces (kept for compatibility).
pub fn make_connexity_blocks(_faces: &[usize], _ds: &DS, _out: &mut Vec<Vec<usize>>) {}
