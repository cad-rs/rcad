use super::*;

impl<'a> super::PaveFiller<'a> {
    pub(crate) fn intersect_cone_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cone1: &ConicalSurface,
        cone2: &ConicalSurface,
    ) {
        use inttools::cone_cone::{ConeConeResult, intersect_cone_cone};
        use inttools::pcurve_derive::{circle_pcurve_on_cone, polyline_pcurve_by_projection};
        use std::f64::consts::TAU;

        // Determine which face is cone1 (for pcurve_on_a/b ordering).
        let cone1_is_f1 = {
            if let Surface3::Cone(c) = &self.ds.faces[f1].surface {
                (c.apex - cone1.apex).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
                    && (c.axis - cone1.axis).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
            } else {
                false
            }
        };

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cone1_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        match intersect_cone_cone(cone1, cone2) {
            ConeConeResult::NoIntersection => (),
            ConeConeResult::Coaxial => {
                self.ds.interferences.push(Interference::FaceFace {
                    f1, f2, curves: vec![], points: vec![],
                });
            }

            ConeConeResult::CoaxialPoint(_pt) => {
                // Single shared apex �?a point contact, not a curve.
            }

            ConeConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
            }

            ConeConeResult::SkewQuartic(branches) => {
                let s1 = Surface3::Cone(*cone1);
                let s2 = Surface3::Cone(*cone2);
                let mut curve_indices = Vec::new();
                for branch in branches {
                    if branch.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(branch[0]);
                    let v_end = self.ds.add_vertex(branch[branch.len() - 1]);
                    let dir = (branch[branch.len() - 1] - branch[0])
                        .normalize_or_zero();
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: branch[0],
                            direction: if dir.length_squared() > 0.5 {
                                dir
                            } else {
                                DVec3::X
                            },
                        }),
                        polyline: branch.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, 1.0],
                        pcurve_on_a: polyline_pcurve_by_projection(&branch, &s1),
                        pcurve_on_b: polyline_pcurve_by_projection(&branch, &s2),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
            curve_extra: crate::bopds::ds::CurveExtra::default(),
        });
                    curve_indices.push(ci);
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
                if !curve_indices.is_empty() {
                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: curve_indices,
                        points: vec![],
                    });
                }
            }

            ConeConeResult::CoaxialCircle(circle) => {
                let pca_cone1 = circle_pcurve_on_cone(&circle, cone1);
                let pcb_cone2 = circle_pcurve_on_cone(&circle, cone2);
                let (pca, pcb) = make_pcurves(pca_cone1, pcb_cone2);

                let pts = sample_circle_arc(&circle, 0.0, TAU, 32);
                let v_start = self.ds.add_vertex(pts[0]);
                let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                let ci = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, TAU],
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
            curve_extra: crate::bopds::ds::CurveExtra::default(),
        });
                self.ds.faces[f1].face_info.curves_sc.insert(ci);
                self.ds.faces[f2].face_info.curves_sc.insert(ci);
                self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![ci],
                    points: vec![],
                });
            }
        }
    }

    // 鈹€鈹€ Torus intersection helpers 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

}
