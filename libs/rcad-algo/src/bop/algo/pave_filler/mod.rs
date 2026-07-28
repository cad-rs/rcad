// OCCT BOPAlgo_PaveFiller — intersection engine.

use crate::bop::algo::{Alert, GlueEnum, Report};
use crate::bop::ds::{
    DS, InterferenceVV, InterferenceVE, InterferenceEE,
    InterferenceVF, InterferenceEF, InterferenceFF, BOPDS_Iterator,
};
use crate::bop::int_tools;
use crate::tolerance::*;
use rcad_kernel::CurveEval;

pub struct PaveFiller<'a> {
    ds: &'a mut DS,
    my_report: Report,
    my_glue: GlueEnum,
    my_fuzzy_value: f64,
    my_run_parallel: bool,
}

impl<'a> PaveFiller<'a> {
    pub fn new(ds: &'a mut DS) -> Self {
        PaveFiller {
            ds,
            my_report: Report::new(),
            my_glue: GlueEnum::GlueOff,
            my_fuzzy_value: 0.0,
            my_run_parallel: false,
        }
    }
    pub fn set_glue(&mut self, enable: bool, tolerance: f64) {
        self.my_glue = if enable { GlueEnum::GlueFull } else { GlueEnum::GlueOff };
        self.my_fuzzy_value = tolerance;
    }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.my_fuzzy_value = v; }
    pub fn set_run_parallel(&mut self, v: bool) { self.my_run_parallel = v; }
    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }
    pub fn report(&self) -> &Report { &self.my_report }

    pub fn perform(&mut self) {
        self.prepare();
        // VV
        self.perform_vv();
        // VE
        self.perform_ve();
        self.update_sd();
        // EE
        self.perform_ee();
        self.update_sd();
        // VF
        self.perform_vf();
        // EF
        self.perform_ef();
        self.update_sd();
        // FF
        self.perform_ff();
    }

    fn prepare(&mut self) {}

    fn update_sd(&self) {
        // Update SD vertices in pave blocks
    }

    // ==================================================================
    // VV
    // ==================================================================
    fn perform_vv(&mut self) {
        let n = self.ds.nb_shapes();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != rcad_kernel::topods::ShapeType::Vertex { continue; }
            let p1 = self.ds.vertex_point_by_idx(i);
            let t1 = self.ds.vertex_tolerance_by_idx(i);
            for j in (i + 1)..n {
                if self.ds.shapes[j].shape_type != rcad_kernel::topods::ShapeType::Vertex { continue; }
                let p2 = self.ds.vertex_point_by_idx(j);
                let t2 = self.ds.vertex_tolerance_by_idx(j);
                let dist = (p1 - p2).length();
                let tol = t1.max(t2);
                if dist <= tol && dist <= 1e-7 {
                    let merged = if t1 >= t2 { i } else { j };
                    self.ds.interf_vv.push(InterferenceVV { v1: i, v2: j, merged_vertex: merged });
                    self.ds.add_interf(i, j);
                }
            }
        }
    }

    // ==================================================================
    // VE
    // ==================================================================
    fn perform_ve(&mut self) {
        for i in 0..self.ds.nb_shapes() {
            if self.ds.shapes[i].shape_type != rcad_kernel::topods::ShapeType::Vertex { continue; }
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            for j in 0..self.ds.nb_shapes() {
                if self.ds.shapes[j].shape_type != rcad_kernel::topods::ShapeType::Edge { continue; }
                let Some(curve) = self.ds.edge_curve(j) else { continue; };
                let (param, proj) = crate::bop::closest_point_on_curve(&curve, pt);
                let dist = (proj - pt).length();
                if dist <= v_tol + 1e-7 {
                    self.ds.interf_ve.push(InterferenceVE {
                        vertex: i, edge: j, param, index_new: 0,
                    });
                    self.ds.add_interf(i, j);
                }
            }
        }
    }

    // ==================================================================
    // EE
    // ==================================================================
    fn perform_ee(&mut self) {
        for i in 0..self.ds.nb_shapes() {
            if self.ds.shapes[i].shape_type != rcad_kernel::topods::ShapeType::Edge { continue; }
            let Some(curve1) = self.ds.edge_curve(i) else { continue; };
            let r1 = self.ds.edge_range(i);
            for j in (i + 1)..self.ds.nb_shapes() {
                if self.ds.shapes[j].shape_type != rcad_kernel::topods::ShapeType::Edge { continue; }
                let Some(curve2) = self.ds.edge_curve(j) else { continue; };
                let r2 = self.ds.edge_range(j);
                let mut ee = int_tools::edge_edge::EdgeEdgeIntersector::new();
                ee.set_edges(i, curve1.clone(), r1, 1e-7, j, curve2.clone(), r2, 1e-7);
                ee.set_fuzzy_value(1e-7);
                ee.perform();
                if !ee.is_done() || ee.common_parts().is_empty() { continue; }
                for cp in ee.common_parts() {
                    let (p1, p2) = if cp.is_edge_type {
                        let r2p = cp.ranges2.first().copied().unwrap_or([0.0, 0.0]);
                        self.ds.interf_ee.push(InterferenceEE {
                            e1: i, e2: j,
                            point: cp.bounding_point1,
                            param1: cp.range1[0],
                            param2: r2p[0],
                            new_vertex: usize::MAX,
                            range1: cp.range1,
                            range2: r2p,
                        });
                        (cp.bounding_point1, cp.bounding_point2)
                    } else {
                        self.ds.interf_ee.push(InterferenceEE {
                            e1: i, e2: j,
                            point: cp.bounding_point1,
                            param1: cp.vertex_param1,
                            param2: cp.vertex_param2,
                            new_vertex: usize::MAX,
                            range1: cp.range1,
                            range2: cp.ranges2.first().copied().unwrap_or([0.0, 0.0]),
                        });
                        (cp.bounding_point1, cp.bounding_point2)
                    };
                    let _ = (p1, p2);
                    self.ds.add_interf(i, j);
                }
            }
        }
    }

    // ==================================================================
    // VF
    // ==================================================================
    fn perform_vf(&mut self) {
        for i in 0..self.ds.nb_shapes() {
            if self.ds.shapes[i].shape_type != rcad_kernel::topods::ShapeType::Vertex { continue; }
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            for j in 0..self.ds.nb_shapes() {
                if self.ds.shapes[j].shape_type != rcad_kernel::topods::ShapeType::Face { continue; }
                let Some(surf) = self.ds.face_surface(j) else { continue; };
                let (uv, proj) = crate::bop::closest_point_on_surface(&surf, pt);
                let dist = (proj - pt).length();
                if dist <= v_tol + 1e-7 {
                    self.ds.interf_vf.push(InterferenceVF {
                        vertex: i, face: j, u: uv.x, v: uv.y,
                        index_new: None,
                    });
                    self.ds.add_interf(i, j);
                }
            }
        }
    }

    // ==================================================================
    // EF
    // ==================================================================
    fn perform_ef(&mut self) {
        for i in 0..self.ds.nb_shapes() {
            if self.ds.shapes[i].shape_type != rcad_kernel::topods::ShapeType::Edge { continue; }
            let Some(curve) = self.ds.edge_curve(i) else { continue; };
            let range = self.ds.edge_range(i);
            for j in 0..self.ds.nb_shapes() {
                if self.ds.shapes[j].shape_type != rcad_kernel::topods::ShapeType::Face { continue; }
                let Some(surf) = self.ds.face_surface(j) else { continue; };
                // Simplified EF check: project curve midpoint to face
                let mid_t = (range[0] + range[1]) * 0.5;
                let mid_pt = curve.point_at(mid_t);
                let (uv, proj) = crate::bop::closest_point_on_surface(&surf, mid_pt);
                let dist = (proj - mid_pt).length();
                if dist <= 1e-6 {
                    self.ds.interf_ef.push(InterferenceEF {
                        edge: i, face: j, point: mid_pt, edge_param: mid_t,
                        new_vertex: usize::MAX,
                    });
                    self.ds.add_interf(i, j);
                }
                let _ = uv;
            }
        }
    }

    // ==================================================================
    // FF — face-face intersection
    // ==================================================================
    fn perform_ff(&mut self) {
        for i in 0..self.ds.nb_shapes() {
            if self.ds.shapes[i].shape_type != rcad_kernel::topods::ShapeType::Face { continue; }
            let Some(surf1) = self.ds.face_surface(i) else { continue; };
            for j in (i + 1)..self.ds.nb_shapes() {
                if self.ds.shapes[j].shape_type != rcad_kernel::topods::ShapeType::Face { continue; }
                let Some(surf2) = self.ds.face_surface(j) else { continue; };
                let ic_curves = int_tools::face_face::intersect_faces(
                    &surf1, &surf2, 1e-7, 1e-7,
                );
                if ic_curves.is_empty() { continue; }
                self.ds.interf_ff.push(InterferenceFF {
                    f1: i, f2: j,
                    curves: ic_curves.into_iter().map(|c| c.tolerance as usize).collect(),
                    points: Vec::new(),
                    tangent_faces: false,
                });
                self.ds.add_interf(i, j);
            }
        }
    }
}
