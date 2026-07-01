use super::*;

impl<'a> super::PaveFiller<'a> {
    pub(crate) fn intersect_cylinder_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cyl: &CylindricalSurface,
        cone: &ConicalSurface,
    ) {
        use inttools::cylinder_cone::{CylinderConeResult, intersect_cylinder_cone};
        use inttools::pcurve_derive::{
            circle_pcurve_on_cone, circle_pcurve_on_cylinder,
            polyline_pcurve_by_projection,
        };
        use std::f64::consts::TAU;

        // Determine which face carries the cylinder (for pcurve_on_a/b ordering).
        let cyl_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Cylinder(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cyl_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        match intersect_cylinder_cone(cyl, cone) {
            CylinderConeResult::NoIntersection => (),

            CylinderConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
            }

            CylinderConeResult::SkewQuartic(branches) => {
                let s1 = Surface3::Cylinder(*cyl);
                let s2 = Surface3::Cone(*cone);
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

            CylinderConeResult::ParallelOffsetPolyline(branches) => {
                let s1 = Surface3::Cylinder(*cyl);
                let s2 = Surface3::Cone(*cone);
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

            CylinderConeResult::CoaxialCircle(circle) => {
                let pca_cyl = circle_pcurve_on_cylinder(&circle, cyl);
                let pcb_cone = circle_pcurve_on_cone(&circle, cone);
                let (pca, pcb) = make_pcurves(pca_cyl, pcb_cone);

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

            CylinderConeResult::CoaxialTwoCircles(c1, c2) => {
                let mut curve_indices = Vec::new();
                for circle in [c1, c2] {
                    let pca_cyl = circle_pcurve_on_cylinder(&circle, cyl);
                    let pcb_cone = circle_pcurve_on_cone(&circle, cone);
                    let (pca, pcb) = make_pcurves(pca_cyl, pcb_cone);

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
                    curve_indices.push(ci);
                }
                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: curve_indices,
                    points: vec![],
                });
            }
        }
    }

    // 鈹€鈹€ Cone �?Cone analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    pub(crate) fn cylinder_face_v_range(&self, face_idx: usize, cyl: &CylindricalSurface) -> [f64; 2] {
        let axis = cyl.axis.normalize();
        let mut v_min = f64::INFINITY;
        let mut v_max = f64::NEG_INFINITY;
        for &ei in &self.ds.faces[face_idx].boundary_edges {
            if let Some(edge) = self.ds.edges.get(ei) {
                if let Some(v) = self.ds.vertices.get(edge.start_vertex) {
                    let proj = (v.point - cyl.origin).dot(axis);
                    v_min = v_min.min(proj);
                    v_max = v_max.max(proj);
                }
                if let Some(v) = self.ds.vertices.get(edge.end_vertex) {
                    let proj = (v.point - cyl.origin).dot(axis);
                    v_min = v_min.min(proj);
                    v_max = v_max.max(proj);
                }
            }
        }
        let r = if v_min.is_infinite() {
            // Fallback: a generous default range
            [-20.0, 20.0]
        } else {
            [v_min, v_max]
        };
        r
    }

    // 鈹€鈹€ Plane �?Cylinder analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    pub(crate) fn intersect_cylinder_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cyl1: &CylindricalSurface,
        cyl2: &CylindricalSurface,
    ) {
        use inttools::cylinder_cylinder::{CylinderCylinderResult, intersect_cylinder_cylinder};
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, ellipse_pcurve_on_cylinder, line_pcurve_on_cylinder, polyline_pcurve_by_projection,
        };
        use std::f64::consts::TAU;

        // Determine which face is cyl1 (for pcurve_on_a/b ordering)
        let cyl1_is_f1 = {
            if let Surface3::Cylinder(c) = &self.ds.faces[f1].surface {
                (c.origin - cyl1.origin).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
                    && (c.axis - cyl1.axis).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
            } else {
                false
            }
        };

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cyl1_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        // Helper: push a circle intersection curve and register it with both faces.
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
                let ci = ds.intersection_curves.len();
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
                ds.faces[f1].face_info.curves_sc.insert(ci);
                ds.faces[f2].face_info.curves_sc.insert(ci);
                ds.faces[f1].face_info.vertices_in.insert(v_start);
                ds.faces[f1].face_info.vertices_in.insert(v_end);
                ds.faces[f2].face_info.vertices_in.insert(v_start);
                ds.faces[f2].face_info.vertices_in.insert(v_end);
                ci
            };

        // Helper: push a line generator intersection and register it.
        let add_line = |ds: &mut DS,
                        line: &Line3,
                        t_range: [f64; 2],
                        pcurve_on_a: Option<Curve2d>,
                        pcurve_on_b: Option<Curve2d>,
                        f1: usize,
                        f2: usize|
         -> usize {
            use rcad_kernel::CurveEval;
            let v_start = ds.add_vertex(Curve3::Line(*line).point_at(t_range[0]));
            let v_end = ds.add_vertex(Curve3::Line(*line).point_at(t_range[1]));
            let ci = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve: Curve3::Line(*line),
                polyline: vec![],
                start_vertex: v_start,
                end_vertex: v_end,
                t_range,
                pcurve_on_a,
                pcurve_on_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
            curve_extra: crate::bopds::ds::CurveExtra::default(),
        });
            ds.faces[f1].face_info.curves_sc.insert(ci);
            ds.faces[f2].face_info.curves_sc.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            ci
        };

        // Helper: push an ellipse intersection and register it.
        let add_ellipse = |ds: &mut DS,
                           ellipse: &Ellipse3,
                           pcurve_on_a: Option<Curve2d>,
                           pcurve_on_b: Option<Curve2d>,
                           f1: usize,
                           f2: usize|
         -> usize {
            let pts = sample_circle_arc(
                &Circle3 {
                    center: ellipse.center,
                    normal: ellipse.normal,
                    radius: ellipse.major_radius.max(ellipse.minor_radius),
                },
                0.0,
                TAU,
                32,
            );
            let v_start = ds.add_vertex(pts[0]);
            let v_end = ds.add_vertex(pts[pts.len() - 1]);
            let ci = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve: Curve3::Ellipse(*ellipse),
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
            ds.faces[f1].face_info.curves_sc.insert(ci);
            ds.faces[f2].face_info.curves_sc.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            ci
        };

        let extent = 20.0_f64;
        let mut curve_indices = Vec::new();

        match intersect_cylinder_cylinder(cyl1, cyl2) {
            CylinderCylinderResult::NoIntersection => return,
            CylinderCylinderResult::Coaxial => {
                // Same-domain coaxial cylinders: record empty-curves FaceFace so
                // the Builder treats this pair as coincident (no intersection to split).
                self.ds.interferences.push(Interference::FaceFace {
                    f1, f2, curves: vec![], points: vec![],
                });
                return;
            }

            CylinderCylinderResult::PerpendicularOffsetCurves {
                cyl1: off_cyl1,
                cyl2: off_cyl2,
                ..
            } => {
                // Perpendicular cylinders with offset (non-intersecting) axes.
                // Parametrization on cyl1's surface:
                //   P(�? = O1 + v(�?*a1 + R1*(cos(�?*U1 + sin(�?*V1)
                //   v(�? = dz �?�?R2�?- (R1路cos(�? - dx)�?
                //
                // Two closed-loop intersection curves per face, one per �?interval:
                //   Loop 1 (胃鈭圼t_low, t_high]): forward branch+  back branch-
                //   Loop 2 (胃鈭圼蟿-t_high, �?t_low]): forward branch+  back branch-
                // Each loop is a single IntersectionCurve whose start/end vertex
                // coincide (same 3D tangent point) �?the boolean builder sees a
                // single closed trim boundary per loop.
                let a1 = off_cyl1.axis.normalize();
                let a2 = off_cyl2.axis.normalize();
                let r1 = off_cyl1.radius;
                let r2 = off_cyl2.radius;
                let r2_sq = r2 * r2;
                let w = off_cyl1.origin - off_cyl2.origin;
                let denom = 1.0 - a1.dot(a2) * a1.dot(a2);
                if denom.abs() < 1e-12 { return; }
                let d1 = a1.dot(w); let d2 = a2.dot(w);
                let t = (a1.dot(a2) * d2 - d1) / denom;
                let s = (d2 - a1.dot(a2) * d1) / denom;
                let conn = (off_cyl1.origin + a1 * t) - (off_cyl2.origin + a2 * s);
                let conn_len = conn.length();
                let u1 = if conn_len < 1e-12 {
                    // Axes intersect (zero or near-zero offset along the
                    // connecting vector).  Pick any direction perpendicular to
                    // a1 �?a1 �?a2 works since the axes are perpendicular.
                    a1.cross(a2).normalize()
                } else {
                    conn / conn_len
                };
                let v1 = a1.cross(u1).normalize();
                let delta = off_cyl2.origin - off_cyl1.origin;
                let dx = delta.dot(u1);
                let dz = delta.dot(a1);
                let cos_min = ((dx - r2) / r1).clamp(-1.0, 1.0);
                let cos_max = ((dx + r2) / r1).clamp(-1.0, 1.0);
                if cos_min > cos_max { return; }
                let t_low = cos_max.acos();
                let t_high = cos_min.acos();

                let surface1 = Surface3::Cylinder(off_cyl1);
                let surface2 = Surface3::Cylinder(off_cyl2);
                let n_per = 9;

                for (t_start, t_end) in [(t_low, t_high), (TAU - t_high, TAU - t_low)] {
                    let n_pts = n_per * 2 + 1; // forward + backward (share the turn-around point)
                    let mut pts: Vec<DVec3> = Vec::with_capacity(n_pts);

                    // Forward: branch = +1, �?= t_start �?t_end
                    for i in 0..=n_per {
                        let theta = t_start + (t_end - t_start) * i as f64 / n_per as f64;
                        let (ct, st) = (theta.cos(), theta.sin());
                        let diff = r1 * ct - dx;
                        let disc = (r2_sq - diff * diff).max(0.0).sqrt();
                        let v_z = dz + disc; // branch sign +1
                        pts.push(off_cyl1.origin + v_z * a1 + r1 * (ct * u1 + st * v1));
                    }
                    // Backward: branch = -1, �?= t_end �?t_start (reversed)
                    for i in 1..=n_per {
                        let theta = t_end - (t_end - t_start) * i as f64 / n_per as f64;
                        let (ct, st) = (theta.cos(), theta.sin());
                        let diff = r1 * ct - dx;
                        let disc = (r2_sq - diff * diff).max(0.0).sqrt();
                        let v_z = dz - disc; // branch sign -1
                        pts.push(off_cyl1.origin + v_z * a1 + r1 * (ct * u1 + st * v1));
                    }
                    if pts.len() < 2 { continue; }

                    let pca = polyline_pcurve_by_projection(&pts, &surface1);
                    let pcb = polyline_pcurve_by_projection(&pts, &surface2);
                    let (pca, pcb) = match (pca, pcb) {
                        (Some(a), Some(b)) => make_pcurves(a, b),
                        _ => continue,
                    };

                    // add_vertex dedup: pts[0] and pts[pts.len()-1] are the same
                    // tangent point (�?t_start, disc=0, both branches coincide).
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: pts[0],
                            direction: (pts[pts.len() - 1] - pts[0]).normalize_or(DVec3::X),
                        }),
                        polyline: pts, start_vertex: v_start, end_vertex: v_end,
                        t_range: [0.0, 1.0], pcurve_on_a: pca, pcurve_on_b: pcb,
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
                    curve_indices.push(ci);
                }
            }

            CylinderCylinderResult::General => {
                // Fall back to numeric marching for skew/oblique axes.
                self.intersect_ff_by_marching(f1, f2);
                return;
            }

            CylinderCylinderResult::OneGeneratorLine(line) => {
                let pca = line_pcurve_on_cylinder(&line, cyl1);
                let pcb = line_pcurve_on_cylinder(&line, cyl2);
                let (pca, pcb) = make_pcurves(pca, pcb);
                let ci = add_line(self.ds, &line, [-extent, extent], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            CylinderCylinderResult::TwoGeneratorLines(l1, l2) => {
                let pca1 = line_pcurve_on_cylinder(&l1, cyl1);
                let pcb1 = line_pcurve_on_cylinder(&l1, cyl2);
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_line(self.ds, &l1, [-extent, extent], pca1, pcb1, f1, f2);

                let pca2 = line_pcurve_on_cylinder(&l2, cyl1);
                let pcb2 = line_pcurve_on_cylinder(&l2, cyl2);
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_line(self.ds, &l2, [-extent, extent], pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::TwoCircles(c1, c2) => {
                // Perpendicular Steinmetz equal-radii: circles in diagonal planes.
                let pca1 = circle_pcurve_on_cylinder(&c1, cyl1);
                let pcb1 = circle_pcurve_on_cylinder(&c1, cyl2);
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_circle(self.ds, &c1, pca1, pcb1, f1, f2);

                let pca2 = circle_pcurve_on_cylinder(&c2, cyl1);
                let pcb2 = circle_pcurve_on_cylinder(&c2, cyl2);
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_circle(self.ds, &c2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::TwoEllipses(e1, e2) => {
                // Perpendicular Steinmetz unequal-radii.
                // Use analytic pcurves on both cylinders (no iterative projection).
                let pca1 = ellipse_pcurve_on_cylinder(&e1, cyl1);
                let pcb1 = ellipse_pcurve_on_cylinder(&e1, cyl2);
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_ellipse(self.ds, &e1, pca1, pcb1, f1, f2);

                let pca2 = ellipse_pcurve_on_cylinder(&e2, cyl1);
                let pcb2 = ellipse_pcurve_on_cylinder(&e2, cyl2);
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_ellipse(self.ds, &e2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::SkewQuartic(branches) => {
                let s1 = Surface3::Cylinder(*cyl1);
                let s2 = Surface3::Cylinder(*cyl2);
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
            }
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
