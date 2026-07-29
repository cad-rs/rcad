// OCCT BRepClass3d_BndBoxTree (BRepClass3d_BndBoxTree.hxx / .cxx)
// BVH tree selectors for point and line queries during solid classification.

use rcad_kernel::base::extrema;
use rcad_kernel::math::bnd::BndBox;
use glam::DVec3;

/// OCCT BRepClass3d_BndBoxTreeSelectorPoint — selects edges/vertices near a point.
pub struct BndBoxTreeSelectorPoint {
    // rcad: simplified — no BVH tree, direct range-based selection
    point: DVec3,
    found: bool,
}

impl BndBoxTreeSelectorPoint {
    pub fn new() -> Self {
        BndBoxTreeSelectorPoint { point: DVec3::ZERO, found: false }
    }

    pub fn set_current_point(&mut self, p: DVec3) { self.point = p; }
    pub fn found(&self) -> bool { self.found }
    pub fn set_found(&mut self, f: bool) { self.found = f; }
    pub fn point(&self) -> DVec3 { self.point }
}

/// OCCT BRepClass3d_BndBoxTreeSelectorLine — selects edges/vertices near a line.
pub struct BndBoxTreeSelectorLine {
    // rcad: simplified — no BVH tree
    line_origin: DVec3,
    line_dir: DVec3,
    max_param: f64,
    // Results
    edge_params: Vec<(usize, f64, f64)>, // (edge_idx, param_on_edge, param_on_line)
    vert_params: Vec<(usize, f64)>,      // (vert_idx, param_on_line)
    is_valid: bool,
}

impl BndBoxTreeSelectorLine {
    pub fn new() -> Self {
        BndBoxTreeSelectorLine {
            line_origin: DVec3::ZERO, line_dir: DVec3::X, max_param: 0.0,
            edge_params: Vec::new(), vert_params: Vec::new(), is_valid: true,
        }
    }

    pub fn set_current_line(&mut self, origin: DVec3, dir: DVec3, max_param: f64) {
        self.line_origin = origin;
        self.line_dir = dir;
        self.max_param = max_param;
    }

    pub fn clear_results(&mut self) {
        self.edge_params.clear();
        self.vert_params.clear();
        self.is_valid = true;
    }

    pub fn is_valid(&self) -> bool { self.is_valid }
    pub fn set_invalid(&mut self) { self.is_valid = false; }

    /// Return edge params collected during selection.
    pub fn edge_params(&self) -> &[(usize, f64, f64)] { &self.edge_params }
    pub fn vert_params(&self) -> &[(usize, f64)] { &self.vert_params }

    pub fn add_edge_param(&mut self, edge_idx: usize, param_on_edge: f64, param_on_line: f64) {
        self.edge_params.push((edge_idx, param_on_edge, param_on_line));
    }
    pub fn add_vert_param(&mut self, vert_idx: usize, param_on_line: f64) {
        self.vert_params.push((vert_idx, param_on_line));
    }

    pub fn nb_edge_params(&self) -> usize { self.edge_params.len() }
    pub fn nb_vert_params(&self) -> usize { self.vert_params.len() }
}
