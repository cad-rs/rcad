use super::*;

impl<'a> super::PaveFiller<'a> {
    pub(crate) fn intersect_sphere_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sphere: &SphericalSurface,
        cone: &ConicalSurface,
    ) {
        use inttools::sphere_cone::{SphereConeResult, intersect_sphere_cone};
        use inttools::pcurve_derive::{
            circle_pcurve_on_cone, fallback_pcurve_by_projection,
            polyline_pcurve_by_projection,
        };
        use std::f64::consts::TAU;

        let sphere_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Sphere(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if sphere_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        let s1 = Surface3::Sphere(*sphere);
        let s2 = Surface3::Cone(*cone);

        match intersect_sphere_cone(sphere, cone) {
            SphereConeResult::NoIntersection => (),

            SphereConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
            }

            SphereConeResult::SingleCircle(circ) => {
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &s1,
                );
                let pcb = circle_pcurve_on_cone(&circ, cone);
                let (pca, pcb) = make_pcurves(pca, pcb);
                let pts = sample_circle_arc(&circ, 0.0, TAU, 32);
                let v_start = self.ds.add_vertex(pts[0]);
                let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                let ci = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circ),
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
                    f1, f2, curves: vec![ci], points: vec![],
                });
            }

            SphereConeResult::TwoCircles(c1, c2) => {
                for circ in [c1, c2] {
                    let pca = fallback_pcurve_by_projection(
                        &Curve3::Circle(circ),
                        &[0.0, TAU],
                        &s1,
                    );
                    let pcb = circle_pcurve_on_cone(&circ, cone);
                    let (pca, pcb) = make_pcurves(pca, pcb);
                    let pts = sample_circle_arc(&circ, 0.0, TAU, 32);
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Circle(circ),
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
                        f1, f2, curves: vec![ci], points: vec![],
                    });
                }
            }

            SphereConeResult::TangentPoint(pt) => {
                let v = self.ds.add_vertex(pt);
                self.ds.faces[f1].face_info.vertices_in.insert(v);
                self.ds.faces[f2].face_info.vertices_in.insert(v);
                self.ds.interferences.push(Interference::FaceFace {
                    f1, f2, curves: vec![], points: vec![v],
                });
            }

            SphereConeResult::Polyline(branches) => {
                let mut curve_indices = Vec::new();
                for branch in branches {
                    if branch.len() < 2 { continue; }
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
                        f1, f2, curves: curve_indices, points: vec![],
                    });
                }
            }
        }
    }
    pub(crate) fn intersect_sphere_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sphere: &SphericalSurface,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, circle_pcurve_on_sphere, polyline_pcurve_by_projection,
        };
        use inttools::sphere_cylinder::{SphereCylinderResult, intersect_sphere_cylinder};
        use std::f64::consts::TAU;

        // Determine which face is the sphere face (for pcurve_on_a/b ordering)
        let sphere_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Sphere(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if sphere_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        // Helper: add one intersection circle to the DS and return its index.
        let add_circle =
            |ds: &mut DS,
             circle: &Circle3,
             pcurve_on_a: Option<Curve2d>,
             pcurve_on_b: Option<Curve2d>,
             f1: usize,
             f2: usize|
             -> usize {
                let pts = sample_circle_arc(circle, 0.0, TAU, 32);
                let v_start = ds.add_vertex(pts[0]);
                let v_end = ds.add_vertex(pts[pts.len() - 1]);
                let curve_idx = ds.intersection_curves.len();
                ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(*circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, TAU],
                    pcurve_on_a,
                    pcurve_on_b,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
            curve_extra: crate::bopds::ds::CurveExtra::default(),
        });
                ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                ds.faces[f1].face_info.vertices_in.insert(v_start);
                ds.faces[f1].face_info.vertices_in.insert(v_end);
                ds.faces[f2].face_info.vertices_in.insert(v_start);
                ds.faces[f2].face_info.vertices_in.insert(v_end);
                curve_idx
            };

        // Closure to compute pcurves for one intersection circle.
        // The intersection circle is always a latitude line on the sphere
        // (�?= acos((h �?h_c) / R)), so `circle_pcurve_on_sphere` is exact
        // here regardless of whether the sphere and cylinder axes are parallel.
        let make_circle_pcurves = |circle: &Circle3| -> (Option<Curve2d>, Option<Curve2d>) {
            let pcurve_sph = circle_pcurve_on_sphere(circle, sphere);
            let pcurve_cyl = circle_pcurve_on_cylinder(circle, cyl);
            make_pcurves(pcurve_sph, pcurve_cyl)
        };

        match intersect_sphere_cylinder(sphere, cyl) {
            SphereCylinderResult::NoIntersection => (),
            SphereCylinderResult::General => {
                // Fall back to numeric marching for the quartic case.
                self.intersect_ff_by_marching(f1, f2);
            }
            SphereCylinderResult::TangentCircle(circle) => {
                let (pca, pcb) = make_circle_pcurves(&circle);
                let ci = add_circle(self.ds, &circle, pca, pcb, f1, f2);
                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![ci],
                    points: vec![],
                });
            }
            SphereCylinderResult::TwoCircles(c1, c2) => {
                let (pca1, pcb1) = make_circle_pcurves(&c1);
                let ci1 = add_circle(self.ds, &c1, pca1, pcb1, f1, f2);
                let (pca2, pcb2) = make_circle_pcurves(&c2);
                let ci2 = add_circle(self.ds, &c2, pca2, pcb2, f1, f2);
                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![ci1, ci2],
                    points: vec![],
                });
            }

            SphereCylinderResult::SkewQuartic(branches) => {
                let s1 = Surface3::Sphere(*sphere);
                let s2 = Surface3::Cylinder(*cyl);
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
        }
    }

    // 鈹€鈹€ Cylinder �?Cylinder analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    pub(crate) fn intersect_sphere_sphere_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sph1: &SphericalSurface,
        sph2: &SphericalSurface,
    ) {
        use inttools::pcurve_derive::fallback_pcurve_by_projection;
        use std::f64::consts::TAU;

        let d_vec = sph2.center - sph1.center;
        let d = d_vec.length();

        // No intersection if disjoint or one contains the other
        if d < TOLERANCE_FLOAT_LOOSE {
            // Concentric spheres: same-domain (same center). Record empty FaceFace.
            self.ds.interferences.push(Interference::FaceFace {
                f1, f2, curves: vec![], points: vec![],
            });
            return;
        }
        if d >= sph1.radius + sph2.radius || d <= (sph1.radius - sph2.radius).abs() {
            return;
        }

        // Distance from sph1 center to the radical plane
        let h = (d * d + sph1.radius * sph1.radius - sph2.radius * sph2.radius) / (2.0 * d);
        let r_circ_sq = sph1.radius * sph1.radius - h * h;
        if r_circ_sq <= 0.0 {
            return; // Tangent or near-tangent
        }
        let r_circ = r_circ_sq.sqrt();

        // Normal of the intersection circle (axis of the radical plane)
        let normal = d_vec.normalize();
        // Center of the intersection circle
        let center = sph1.center + normal * h;

        let circle = Circle3 {
            center,
            normal,
            radius: r_circ,
        };

        let curve3 = Curve3::Circle(circle);
        let t_range = [0.0_f64, TAU];
        let pcurve_a = fallback_pcurve_by_projection(&curve3, &t_range, &Surface3::Sphere(*sph1));
        let pcurve_b = fallback_pcurve_by_projection(&curve3, &t_range, &Surface3::Sphere(*sph2));

        let pts = sample_circle_arc(&circle, 0.0, TAU, 32);
        if pts.len() < 2 {
            return;
        }

        let v_start = self.ds.add_vertex(pts[0]);
        let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

        let curve_idx = self.ds.intersection_curves.len();
        self.ds.intersection_curves.push(IntersectionCurve {
            curve: curve3,
            polyline: vec![],
            start_vertex: v_start,
            end_vertex: v_end,
            t_range: [0.0, TAU],
            pcurve_on_a: Some(pcurve_a),
            pcurve_on_b: Some(pcurve_b),
            geom_tol: crate::tolerance::TOLERANCE_ABS,
        pave_blocks: Vec::new(),
            curve_extra: crate::bopds::ds::CurveExtra::default(),
        });

        self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
        self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
        self.ds.faces[f1].face_info.vertices_in.insert(v_start);
        self.ds.faces[f1].face_info.vertices_in.insert(v_end);
        self.ds.faces[f2].face_info.vertices_in.insert(v_start);
        self.ds.faces[f2].face_info.vertices_in.insert(v_end);

        self.ds.interferences.push(Interference::FaceFace {
            f1,
            f2,
            curves: vec![curve_idx],
            points: vec![],
        });
    }

    // 鈹€鈹€ Sphere �?Cylinder analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

}
