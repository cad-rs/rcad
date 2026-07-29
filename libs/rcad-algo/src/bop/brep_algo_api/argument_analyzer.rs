// OCCT BOPAlgo_ArgumentAnalyzer — input validation for Boolean Operations.
//
// OCCT BOPAlgo_ArgumentAnalyzer.cxx L1-1015 / .hxx.

use crate::bop::algo::checker_si::CheckerSI;
use crate::bop::ds::DS;
use rcad_kernel::topods::{self, Orientation, TShape, ShapeType};
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
        self.prepare();
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

    fn push(&mut self, s: CheckStatus) { self.my_result.push(CheckResult { check_status: s }); }
    fn should_stop(&self) -> bool { self.my_stop_on_first && !self.my_result.is_empty() }
    fn shape_opt(&self, i: usize) -> Option<Shape> {
        if i == usize::MAX { None } else { Some(Shape::synthetic(i, Orientation::Forward)) }
    }

    // ── Prepare L115-126 ────────────────────────────────────────────
    fn prepare(&mut self) {
        // OCCT L117-125: BOPTools_AlgoTools3D::IsEmptyShape
        // rcad: DS shape emptiness not yet queried.
    }

    // ── TestTypes L275-352 ──────────────────────────────────────────
    fn test_types(&mut self) {
        let is_s1 = self.my_shape1 == usize::MAX;
        let is_s2 = self.my_shape2 == usize::MAX;
        if is_s1 && is_s2 { self.push(CheckStatus::BadType); return; }
        if (is_s1 && !is_s2) || (!is_s1 && is_s2) {
            let empty = if is_s1 { self.my_empty2 } else { self.my_empty1 };
            if empty || self.my_operation == 0 { self.push(CheckStatus::BadType); }
            return;
        }
        if self.my_empty1 || self.my_empty2 { self.push(CheckStatus::BadType); return; }
        // OCCT L330-350: operation-specific dimension checks
        // BOPTools_AlgoTools::Dimensions(minDim, maxDim) — rcad: use DS shape exploration
        if self.my_operation != 0 && self.my_operation != 1 { // not UNKNOWN or COMMON
            // rcad: dimension checks not yet implemented
        }
    }

    // ── TestSelfInterferences L356-445 ──────────────────────────────
    fn test_self_interferences(&mut self) {
        for ii in 0..2 {
            let shape = if ii == 0 { self.my_shape1 } else { self.my_shape2 };
            if shape == usize::MAX { continue; }
            let mut ds = DS::new();
            ds.set_arguments(vec![Shape::synthetic(shape, Orientation::Forward)]);
            ds.init(1e-7);
            let mut checker = CheckerSI::new();
            checker.perform(&mut ds);
            // L390-426: iterate interferences, skip new shapes
            let pairs: Vec<_> = ds.interf_tb.iter().copied().collect();
            for (n1, n2) in &pairs {
                if ds.is_new_shape(*n1) || ds.is_new_shape(*n2) { continue; }
                self.push(CheckStatus::SelfIntersect);
                if self.should_stop() { return; }
            }
            // L428-443: if checker had errors, add OperationAborted
        }
    }

    // ── TestSmallEdge L449-567 ──────────────────────────────────────
    fn test_small_edge(&mut self) {
        let a_ctx = crate::bop::int_tools::context::IntToolsContext::new();
        for i in 0..2 {
            let shape_i = if i == 0 { self.my_shape1 } else { self.my_shape2 };
            if shape_i == usize::MAX { continue; }
            // OCCT TopExp_Explorer(aS, TopAbs_EDGE)
            // rcad: iterate DS shapes that are edges
            let mut ds = DS::new();
            ds.set_arguments(vec![Shape::synthetic(shape_i, Orientation::Forward)]);
            ds.init(1e-7);
            let n_src = ds.nb_source_shapes();
            for ei in 0..n_src {
                let si = ds.shape_info(ei);
                if si.shape_type != ShapeType::Edge { continue; }
                // OCCT check degenerated
                let ed_data = match &*si.shape.data {
                    TShape::Edge(e) => e,
                    _ => continue,
                };
                if ed_data.degenerated { continue; }
                let tol_e = ed_data.tolerance;
                // OCCT IsMicroEdge: checks if edge length < tolerance
                let range = ed_data.range;
                let is_micro = (range[1] - range[0]).abs() < tol_e * 2.0;
                if !is_micro { continue; }
                // For SECTION: check if edge vertices lie on other shape
                let mut keep = true;
                if self.my_operation == 5 && self.my_shape2 != usize::MAX {
                    // OCCT L481-537: project edge vertices onto other shape
                    // check distance vs tolerance sum
                    let _ = &a_ctx;
                    keep = false; // rcad stub
                }
                if keep {
                    self.push(CheckStatus::TooSmallEdge);
                    if self.should_stop() { return; }
                }
            }
        }
    }

    // ── TestRebuildFace L571-672 ────────────────────────────────────
    fn test_rebuild_face(&mut self) {
        if self.my_operation == 5 || self.my_operation == 0 { return; }
        for i in 0..2 {
            let shape_i = if i == 0 { self.my_shape1 } else { self.my_shape2 };
            if shape_i == usize::MAX { continue; }
            let mut ds = DS::new();
            ds.set_arguments(vec![Shape::synthetic(shape_i, Orientation::Forward)]);
            ds.init(1e-7);
            let n_src = ds.nb_source_shapes();
            for fi in 0..n_src {
                let si = ds.shape_info(fi);
                if si.shape_type != ShapeType::Face { continue; }
                // OCCT L592-645: count internal edges, run BOPAlgo_BuilderFace,
                // verify 1 area and matching edge count
                // rcad: BuilderFace not yet integrated
                let _ = si;
            }
        }
    }

    // ── TestTangent L676-679 ────────────────────────────────────────
    fn test_tangent(&mut self) { /* OCCT: not implemented */ }

    // ── TestMergeSubShapes L683-878 ─────────────────────────────────
    fn test_merge_sub_shapes(&mut self, the_type: ShapeType) {
        if self.my_shape1 == usize::MAX || self.my_shape2 == usize::MAX { return; }
        if self.my_empty1 || self.my_empty2 { return; }
        let status = match the_type {
            ShapeType::Vertex => CheckStatus::IncompatibilityOfVertex,
            ShapeType::Edge => CheckStatus::IncompatibilityOfEdge,
            ShapeType::Face => CheckStatus::IncompatibilityOfFace,
            _ => return,
        };
        // OCCT L714-740: collect sub-shapes from both shapes
        let mut ds = DS::new();
        ds.set_arguments(vec![
            Shape::synthetic(self.my_shape1, Orientation::Forward),
            Shape::synthetic(self.my_shape2, Orientation::Forward),
        ]);
        ds.init(1e-7);
        let n_src = ds.nb_source_shapes();
        let mut seq1: Vec<usize> = Vec::new();
        let mut seq2: Vec<usize> = Vec::new();
        // First shape sub-shapes
        for i in 0..n_src {
            if ds.shape_info(i).shape_type == the_type && ds.rank(i) == 0 {
                seq1.push(i);
            }
        }
        // Second shape sub-shapes
        for i in 0..n_src {
            if ds.shape_info(i).shape_type == the_type && ds.rank(i) == 1 {
                seq2.push(i);
            }
        }
        // OCCT L752-835: compare shape1 sub-shapes with shape2
        for &s1 in &seq1 {
            let mut matches: Vec<usize> = Vec::new();
            for &s2 in &seq2 {
                let eq = match the_type {
                    ShapeType::Vertex => {
                        // distance vs tolerance sum
                        let p1 = ds.vertex_point_by_idx(s1);
                        let p2 = ds.vertex_point_by_idx(s2);
                        let t1 = ds.vertex_tolerance_by_idx(s1);
                        let t2 = ds.vertex_tolerance_by_idx(s2);
                        (p1 - p2).length() <= t1 + t2
                    }
                    ShapeType::Edge => {
                        // OCCT IntTools_EdgeEdge for coincidence
                        // rcad: simplified — check if bounding vertices are close
                        let si1 = ds.shape_info(s1);
                        let si2 = ds.shape_info(s2);
                        let subs1: Vec<usize> = si1.sub_shapes.iter().copied().collect();
                        let subs2: Vec<usize> = si2.sub_shapes.iter().copied().collect();
                        subs1.len() == subs2.len() && subs1.iter().zip(&subs2).all(|(a, b)| a == b)
                    }
                    _ => false,
                };
                if eq { matches.push(s2); }
            }
            if matches.len() > 1 {
                self.push(status);
                if self.should_stop() { return; }
            }
        }
        // OCCT L838-877: reverse check (shape2 vs shape1)
        for &s2 in &seq2 {
            let mut matches: Vec<usize> = Vec::new();
            for &s1 in &seq1 {
                let eq = match the_type {
                    ShapeType::Vertex => {
                        let p1 = ds.vertex_point_by_idx(s1);
                        let p2 = ds.vertex_point_by_idx(s2);
                        let t1 = ds.vertex_tolerance_by_idx(s1);
                        let t2 = ds.vertex_tolerance_by_idx(s2);
                        (p1 - p2).length() <= t1 + t2
                    }
                    _ => false,
                };
                if eq { matches.push(s1); }
            }
            if matches.len() > 1 {
                self.push(status);
                if self.should_stop() { return; }
            }
        }
    }

    fn test_merge_vertex(&mut self) { self.test_merge_sub_shapes(ShapeType::Vertex); }
    fn test_merge_edge(&mut self) { self.test_merge_sub_shapes(ShapeType::Edge); }

    // ── TestContinuity L896-958 ────────────────────────────────────
    fn test_continuity(&mut self) {
        for i in 0..2 {
            let shape_i = if i == 0 { self.my_shape1 } else { self.my_shape2 };
            if shape_i == usize::MAX { continue; }
            let mut ds = DS::new();
            ds.set_arguments(vec![Shape::synthetic(shape_i, Orientation::Forward)]);
            ds.init(1e-7);
            let n_src = ds.nb_source_shapes();
            for j in 0..n_src {
                let si = ds.shape_info(j);
                // OCCT L912-925: check edges with C0 continuity
                if si.shape_type == ShapeType::Edge {
                    // rcad: curve continuity not yet available on Curve3
                }
                // OCCT L927-936: check faces with C0 continuity
                if si.shape_type == ShapeType::Face {
                    // rcad: surface continuity not yet available on Surface3
                }
            }
        }
    }

    // ── TestCurveOnSurface L962-1015 ────────────────────────────────
    fn test_curve_on_surface(&mut self) {
        for i in 0..2 {
            let shape_i = if i == 0 { self.my_shape1 } else { self.my_shape2 };
            if shape_i == usize::MAX { continue; }
            let mut ds = DS::new();
            ds.set_arguments(vec![Shape::synthetic(shape_i, Orientation::Forward)]);
            ds.init(1e-7);
            let n_src = ds.nb_source_shapes();
            for fi in 0..n_src {
                let fsi = ds.shape_info(fi);
                if fsi.shape_type != ShapeType::Face { continue; }
                // OCCT L978-1013: for each edge on face, call ComputeTolerance
                for &ei in &fsi.sub_shapes {
                    let esi = ds.shape_info(ei);
                    if esi.shape_type != ShapeType::Edge { continue; }
                    // BOPTools_AlgoTools::ComputeTolerance(aF, aE, aD, aT)
                    // rcad: not yet implemented
                    let _ = esi;
                }
            }
        }
    }
}

impl Default for ArgumentAnalyzer {
    fn default() -> Self { Self::new() }
}
