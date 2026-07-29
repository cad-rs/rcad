// OCCT BOPAlgo_ArgumentAnalyzer — input validation for Boolean Operations.
//
// OCCT BOPAlgo_ArgumentAnalyzer.cxx L1-1015 / .hxx.

use crate::bop::algo::checker_si::CheckerSI;
use crate::bop::ds::DS;
use rcad_kernel::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

/// OCCT BOPAlgo_ArgumentAnalyzer — checks shape validity for Boolean Ops.
pub struct ArgumentAnalyzer {
    my_shape1: usize,
    my_shape2: usize,
    my_stop_on_first: bool,
    my_operation: i32,
    my_argument_type_mode: bool,
    my_self_inter_mode: bool,
    my_small_edge_mode: bool,
    my_rebuild_face_mode: bool,
    my_tangent_mode: bool,
    my_merge_vertex_mode: bool,
    my_merge_edge_mode: bool,
    my_continuity_mode: bool,
    my_curve_on_surface_mode: bool,
    my_empty1: bool,
    my_empty2: bool,
    my_result: Vec<CheckResult>,
}

#[derive(Debug, Clone)]
pub struct CheckResult {
    pub check_status: CheckStatus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckStatus {
    BadType,
    SelfIntersect,
    TooSmallEdge,
    NonRecoverableFace,
    IncompatibilityOfVertex,
    IncompatibilityOfEdge,
    IncompatibilityOfFace,
    GeomAbsC0,
    InvalidCurveOnSurface,
    OperationAborted,
    CheckUnknown,
}

impl ArgumentAnalyzer {
    pub fn new() -> Self {
        ArgumentAnalyzer {
            my_shape1: usize::MAX,
            my_shape2: usize::MAX,
            my_stop_on_first: false,
            my_operation: 0,
            my_argument_type_mode: false,
            my_self_inter_mode: false,
            my_small_edge_mode: false,
            my_rebuild_face_mode: false,
            my_tangent_mode: false,
            my_merge_vertex_mode: false,
            my_merge_edge_mode: false,
            my_continuity_mode: false,
            my_curve_on_surface_mode: false,
            my_empty1: false,
            my_empty2: false,
            my_result: Vec::new(),
        }
    }

    pub fn set_shape1(&mut self, s: usize) { self.my_shape1 = s; }
    pub fn set_shape2(&mut self, s: usize) { self.my_shape2 = s; }
    pub fn get_shape1(&self) -> usize { self.my_shape1 }
    pub fn get_shape2(&self) -> usize { self.my_shape2 }
    pub fn operation_type(&mut self) -> &mut i32 { &mut self.my_operation }
    pub fn stop_on_first_faulty(&mut self) -> &mut bool { &mut self.my_stop_on_first }
    pub fn argument_type_mode(&mut self) -> &mut bool { &mut self.my_argument_type_mode }
    pub fn self_inter_mode(&mut self) -> &mut bool { &mut self.my_self_inter_mode }
    pub fn small_edge_mode(&mut self) -> &mut bool { &mut self.my_small_edge_mode }
    pub fn rebuild_face_mode(&mut self) -> &mut bool { &mut self.my_rebuild_face_mode }
    pub fn tangent_mode(&mut self) -> &mut bool { &mut self.my_tangent_mode }
    pub fn merge_vertex_mode(&mut self) -> &mut bool { &mut self.my_merge_vertex_mode }
    pub fn merge_edge_mode(&mut self) -> &mut bool { &mut self.my_merge_edge_mode }
    pub fn continuity_mode(&mut self) -> &mut bool { &mut self.my_continuity_mode }
    pub fn curve_on_surface_mode(&mut self) -> &mut bool { &mut self.my_curve_on_surface_mode }

    /// OCCT Perform() L130-257.
    pub fn perform(&mut self) {
        self.my_result.clear();
        self.prepare();                                              // L142
        if self.my_argument_type_mode       { self.test_types(); }
        if self.my_self_inter_mode          { self.test_self_interferences(); if self.should_stop() { return; } }
        if self.my_small_edge_mode          { self.test_small_edge(); if self.should_stop() { return; } }
        if self.my_rebuild_face_mode        { self.test_rebuild_face(); if self.should_stop() { return; } }
        if self.my_tangent_mode             { self.test_tangent(); if self.should_stop() { return; } }
        if self.my_merge_vertex_mode        { self.test_merge_vertex(); if self.should_stop() { return; } }
        if self.my_merge_edge_mode          { self.test_merge_edge(); if self.should_stop() { return; } }
        if self.my_continuity_mode          { self.test_continuity(); if self.should_stop() { return; } }
        if self.my_curve_on_surface_mode    { self.test_curve_on_surface(); }
    }

    pub fn has_faulty(&self) -> bool { !self.my_result.is_empty() }
    pub fn get_check_result(&self) -> &[CheckResult] { &self.my_result }

    fn prepare(&mut self) {
        // OCCT L115-126: BOPTools_AlgoTools3D::IsEmptyShape
        // rcad: not yet implemented.
    }

    fn test_types(&mut self) {
        // OCCT L275-352: shape type checks against operation
        let is_s1 = self.my_shape1 == usize::MAX;
        let is_s2 = self.my_shape2 == usize::MAX;
        if is_s1 && is_s2 { self.push(CheckStatus::BadType); return; }
        if (is_s1 && !is_s2) || (!is_s1 && is_s2) {
            let empty = if is_s1 { self.my_empty2 } else { self.my_empty1 };
            if empty || self.my_operation == 0 { self.push(CheckStatus::BadType); }
            return;
        }
        if self.my_empty1 || self.my_empty2 { self.push(CheckStatus::BadType); return; }
        // L330: operation-specific dimension checks (requires BOPTools_AlgoTools::Dimensions)
    }

    fn test_self_interferences(&mut self) {
        // OCCT L356-445: run CheckerSI, iterate interferences
        for ii in 0..2 {
            let shape = if ii == 0 { self.my_shape1 } else { self.my_shape2 };
            if shape == usize::MAX { continue; }
            let mut ds = DS::new();
            ds.set_arguments(vec![Shape::synthetic(shape, rcad_kernel::topods::Orientation::Forward)]);
            ds.init(1e-7);
            let mut checker = CheckerSI::new();
            checker.perform(&mut ds);
            // OCCT L390-426: iterate ds.Interferences(), add SelfIntersect results
            let a_it: Vec<_> = ds.interf_tb.iter().copied().collect();
            for (n1, n2) in &a_it {
                if ds.is_new_shape(*n1) || ds.is_new_shape(*n2) { continue; }
                self.push(CheckStatus::SelfIntersect);
                if self.should_stop() { return; }
            }
            // OCCT L428-443: if checker has errors, add OperationAborted
        }
    }

    fn test_small_edge(&mut self) {
        // OCCT L449-567: iterate edges, check IsMicroEdge
        // rcad: requires TopExp_Explorer equivalent and IsMicroEdge
    }

    fn test_rebuild_face(&mut self) {
        if self.my_operation == 5 || self.my_operation == 0 { return; }
        // OCCT L571-672: iterate faces, run BOPAlgo_BuilderFace,
        // verify exactly 1 area with matching edge count
        // rcad: requires face iteration from DS
    }

    fn test_tangent(&mut self) { /* OCCT: not implemented */ }

    fn test_merge_sub_shapes(&mut self, _the_type: u8) {
        // OCCT L683-878: compare sub-shapes between two shapes
        if self.my_shape1 == usize::MAX || self.my_shape2 == usize::MAX { return; }
        if self.my_empty1 || self.my_empty2 { return; }
        // rcad: sub-shape iteration and comparison not yet wired
    }

    fn test_merge_vertex(&mut self) { self.test_merge_sub_shapes(0); }
    fn test_merge_edge(&mut self) { self.test_merge_sub_shapes(1); }

    fn test_continuity(&mut self) {
        // OCCT L896-958: check C0 edges/faces
        // rcad: continuity query not yet on curve/surface types
    }

    fn test_curve_on_surface(&mut self) {
        // OCCT L962-1015: ComputeTolerance for face/edge pairs
        // rcad: not yet implemented
    }

    fn push(&mut self, s: CheckStatus) { self.my_result.push(CheckResult { check_status: s }); }
    fn should_stop(&self) -> bool { self.my_stop_on_first && !self.my_result.is_empty() }
}

impl Default for ArgumentAnalyzer {
    fn default() -> Self { Self::new() }
}
