// OCCT BOPAlgo_PaveFiller 閳?intersection engine.

use crate::bop::algo::{Alert, GlueEnum, Report};
use crate::bop::ds::{
    DS, InterferenceVV, InterferenceVE, InterferenceEE,
    InterferenceVF, InterferenceEF, InterferenceFF, BOPDS_Iterator,
};
use crate::bop::int_tools;

use rcad_kernel::CurveEval;
use rcad_kernel::topods::ShapeType;

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
            ds, my_report: Report::new(),
            my_glue: GlueEnum::GlueOff, my_fuzzy_value: 0.0,
            my_run_parallel: false,
        }
    }
    pub fn set_glue(&mut self, enable: bool, tolerance: f64) {
        self.my_glue = if enable { GlueEnum::GlueFull } else { GlueEnum::GlueOff };
        self.my_fuzzy_value = tolerance;
    }
    pub fn fuzzy_value(&self) -> f64 { self.my_fuzzy_value }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.my_fuzzy_value = v; }
    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }
    pub fn report(&self) -> &Report { &self.my_report }

    pub fn perform(&mut self) {
        self.prepare();
        self.perform_vv();
        self.perform_ve();
        self.perform_ee();
        self.perform_vf();
        self.perform_ef();
        self.perform_ff();
    }

    fn prepare(&mut self) {}

    // 閳光偓閳光偓 VV 閳光偓閳光偓
    fn perform_vv(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_vv: Vec<InterferenceVV> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex { continue; }
            let p1 = self.ds.vertex_point_by_idx(i);
            let t1 = self.ds.vertex_tolerance_by_idx(i);
            for j in (i + 1)..n {
                if self.ds.shapes[j].shape_type != ShapeType::Vertex { continue; }
                let p2 = self.ds.vertex_point_by_idx(j);
                let t2 = self.ds.vertex_tolerance_by_idx(j);
                let dist = (p1 - p2).length();
                let tol = t1.max(t2);
                if dist <= tol.max(1e-7) {
                    let merged = if t1 >= t2 { i } else { j };
                    new_vv.push(InterferenceVV { v1: i, v2: j, merged_vertex: merged });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_vv.extend(new_vv);
    }

    // 閳光偓閳光偓 VE 閳光偓閳光偓
    fn perform_ve(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_ve: Vec<InterferenceVE> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex { continue; }
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            for j in 0..n {
                if self.ds.shapes[j].shape_type != ShapeType::Edge { continue; }
                let Some(curve) = self.ds.edge_curve(j) else { continue; };
                let (param, proj) = crate::bop::closest_point_on_curve(curve, pt);
                let dist = (proj - pt).length();
                if dist <= v_tol + 1e-7 {
                    new_ve.push(InterferenceVE { vertex: i, edge: j, param, index_new: 0 });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_ve.extend(new_ve);
    }

    // 閳光偓閳光偓 EE 閳光偓閳光偓
    fn perform_ee(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_ee: Vec<InterferenceEE> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Edge { continue; }
            let Some(c1) = self.ds.edge_curve(i).cloned() else { continue; };
            let r1 = self.ds.edge_range(i);
            for j in (i + 1)..n {
                if self.ds.shapes[j].shape_type != ShapeType::Edge { continue; }
                let Some(c2) = self.ds.edge_curve(j).cloned() else { continue; };
                let r2 = self.ds.edge_range(j);
                let mut ee = int_tools::edge_edge::EdgeEdgeIntersector::new();
                ee.set_edges(i, r1, j, r2, self.ds);
                ee.set_fuzzy_value(1e-7);
                ee.perform();
                if !ee.is_done() || ee.common_parts().is_empty() { continue; }
                for cp in ee.common_parts() {
                    let r2p = cp.ranges2.first().copied().unwrap_or([0.0, 0.0]);
                    new_ee.push(InterferenceEE {
                        e1: i, e2: j, point: cp.bounding_point1,
                        param1: cp.range1[0], param2: r2p[0],
                        new_vertex: usize::MAX, range1: cp.range1, range2: r2p,
                    });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_ee.extend(new_ee);
    }

    // 閳光偓閳光偓 VF 閳光偓閳光偓
    fn perform_vf(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_vf: Vec<InterferenceVF> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex { continue; }
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            for j in 0..n {
                if self.ds.shapes[j].shape_type != ShapeType::Face { continue; }
                let Some(surf) = self.ds.face_surface(j) else { continue; };
                let (uv, proj) = crate::bop::closest_point_on_surface(&surf, pt);
                let dist = (proj - pt).length();
                if dist <= v_tol + 1e-7 {
                    new_vf.push(InterferenceVF { vertex: i, face: j, u: uv.x, v: uv.y, index_new: None });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_vf.extend(new_vf);
    }

    // 閳光偓閳光偓 EF 閳光偓閳光偓
    fn perform_ef(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_ef: Vec<InterferenceEF> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Edge { continue; }
            let Some(curve) = self.ds.edge_curve(i) else { continue; };
            let range = self.ds.edge_range(i);
            let mid_t = (range[0] + range[1]) * 0.5;
            let mid_pt = curve.point_at(mid_t);
            for j in 0..n {
                if self.ds.shapes[j].shape_type != ShapeType::Face { continue; }
                let Some(surf) = self.ds.face_surface(j) else { continue; };
                let (uv, proj) = crate::bop::closest_point_on_surface(&surf, mid_pt);
                let dist = (proj - mid_pt).length();
                if dist <= 1e-6 {
                    new_ef.push(InterferenceEF { edge: i, face: j, point: mid_pt, edge_param: mid_t, new_vertex: usize::MAX });
                    self.ds.add_interf(i, j);
                }
                let _ = uv;
            }
        }
        self.ds.interf_ef.extend(new_ef);
    }

    // 閳光偓閳光偓 FF 閳光偓閳光偓
    fn perform_ff(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_ff: Vec<InterferenceFF> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Face { continue; }
            let Some(s1) = self.ds.face_surface(i) else { continue; };
            for j in (i + 1)..n {
                if self.ds.shapes[j].shape_type != ShapeType::Face { continue; }
                let Some(s2) = self.ds.face_surface(j) else { continue; };
                let mut ff = int_tools::face_face::FaceFace::new();
                ff.set_surfaces(s1.clone(), s2.clone());
                ff.set_tolerances(1e-7, 1e-7);
                ff.perform();
                if !ff.has_intersection() { continue; }
                new_ff.push(InterferenceFF {
                    f1: i, f2: j,
                    curves: ff.make_curves().iter().map(|c| c.tolerance as usize).collect(),
                    points: Vec::new(),
                    tangent_faces: false,
                });
                self.ds.add_interf(i, j);
            }
        }
        self.ds.interf_ff.extend(new_ff);
    }
}
