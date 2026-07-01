use super::*;

impl<'a> super::PaveFiller<'a> {
    pub(crate) fn intersect_plane_plane_faces(&mut self, f1: usize, f2: usize, p1: &Plane, p2: &Plane) {
        use inttools::pcurve_derive::line_pcurve_on_plane;

        match inttools::plane_plane::intersect_plane_plane(p1, p2) {
            inttools::plane_plane::PlanePlaneResult::Parallel => {
            }
            inttools::plane_plane::PlanePlaneResult::Coincident => {
                self.handle_coplanar_faces(f1, f2, p1);
            }
            inttools::plane_plane::PlanePlaneResult::Line(line) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);
                let clip_tol = self.ff_tol(f1, f2);

                let ranges1 =
                    inttools::edge_face::clip_line_to_polygon_with_tol(&line, p1, &verts1, clip_tol);
                let ranges2 =
                    inttools::edge_face::clip_line_to_polygon_with_tol(&line, p2, &verts2, clip_tol);

                for &(t1_min, t1_max) in &ranges1 {
                    for &(t2_min, t2_max) in &ranges2 {
                        let t_min = t1_min.max(t2_min);
                        let t_max = t1_max.min(t2_max);
                        // Keep strict: overlap length along the intersection line is parametric, not V鈥揤
                        // coincidence �?tying this to `fuzzy_tol` can change sphere鈥揵ox trims and area.
                        if t_max - t_min < TOLERANCE_ABS {
                            continue;
                        }

                        let p_start = line.origin + line.direction * t_min;
                        let p_end = line.origin + line.direction * t_max;

                        let v_start = self.ds.add_vertex(p_start);
                        let v_end = self.ds.add_vertex(p_end);

                        let curve_idx = self.ds.intersection_curves.len();
                        let pca = line_pcurve_on_plane(&line, p1);
                        let pcb = line_pcurve_on_plane(&line, p2);
                        self.ds.intersection_curves.push(IntersectionCurve {
                            curve: Curve3::Line(line),
                            polyline: vec![],
                            start_vertex: v_start,
                            end_vertex: v_end,
                            t_range: [t_min, t_max],
                            pcurve_on_a: Some(pca),
                            pcurve_on_b: Some(pcb),
                            geom_tol: crate::tolerance::TOLERANCE_ABS,
                        pave_blocks: Vec::new(),
            curve_extra: crate::bopds::ds::CurveExtra::default(),
        });

                        self.ds.interferences.push(Interference::FaceFace {
                            f1,
                            f2,
                            curves: vec![curve_idx],
                            points: vec![],
                        });

                        self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                        self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                        self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                        self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                        self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                        self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                    }
                }
            }
        }
    }

    pub(crate) fn handle_coplanar_faces(&mut self, f1: usize, f2: usize, plane: &Plane) {
        let verts1 = self.ds.face_boundary_points(f1);
        let verts2 = self.ds.face_boundary_points(f2);

        let result = inttools::coplanar::analyze_coplanar_faces(&verts1, &verts2, plane);

        if !result.overlap.is_empty() {
            // �?OCCT-aligned: create IC for each overlap edge (BOPAlgo_PaveFiller_6.cxx:285-622)
            let plane1 = self.ds.face_plane(f1);
            let plane2 = self.ds.face_plane(f2);

            for overlap_poly in &result.overlap {
                if overlap_poly.len() < 3 { continue; }
                for i in 0..overlap_poly.len() {
                    let j = (i + 1) % overlap_poly.len();
                    let p_start = overlap_poly[i];
                    let p_end = overlap_poly[j];
                    if (p_end - p_start).length_squared() < TOLERANCE_ABS_SQ { continue; }

                    let v_start = self.ds.add_vertex(p_start);
                    let v_end = self.ds.add_vertex(p_end);
                    let dir = (p_end - p_start).normalize();
                    let len = (p_end - p_start).length();
                    let line = Line3 { origin: p_start, direction: dir };

                    let pca = inttools::pcurve_derive::line_pcurve_on_plane(&line, &plane1);
                    let pcb = inttools::pcurve_derive::line_pcurve_on_plane(&line, &plane2);

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(line),
                        polyline: vec![],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, len],
                        pcurve_on_a: Some(pca),
                        pcurve_on_b: Some(pcb),
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
                }
            }

            // Keep existing same_domain_overlaps for backward compatibility
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: vec![],
                points: vec![],
            });
            if let Some(overlap) = result.overlap.into_iter().max_by_key(|poly| poly.len()) {
                self.ds.same_domain_overlaps.push((f1, f2, overlap));
            }
        }
    }

    // 鈹€鈹€ Plane �?Sphere analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€


    /// OCCT: plane-sphere intersection
    /// OCCT: plane-sphere intersection
    pub(crate) fn intersect_plane_sphere_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        sphere: &SphericalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_plane, circle_pcurve_on_sphere,
        };
        use inttools::plane_sphere::{PlaneSphereResult, intersect_plane_sphere};

        // Determine which face carries the plane (for correct pcurve_on_a/b assignment)
        let plane_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Plane(_));

        let ps_result = intersect_plane_sphere(plane, sphere);
        match ps_result {
            PlaneSphereResult::NoIntersection => {}
            PlaneSphereResult::TangentPoint(pt) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);
                let tff = self.ff_tol(f1, f2);
                if inttools::edge_face::point_in_planar_face_with_tol(pt, plane, &verts1, tff)
                    && point_in_sphere_face(pt, &verts2, self.ds)
                {
                    let v = self.ds.add_vertex(pt);
                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![],
                        points: vec![v],
                    });
                }
            }
            PlaneSphereResult::Circle(circle) => {
                // OCCT IntTools_FaceFace stores the full analytic curve without
                // pre-clipping to face boundaries.  OCCT PaveFiller clips via
                // PutBoundPaveOnCurve during MakeBlocks.
                if circle.radius <= TOLERANCE_MESH_LEGACY + TOLERANCE_ABS { return; }
                let (effective_t0, effective_t1) = (0.0_f64, std::f64::consts::TAU);

                let pcurve_plane = circle_pcurve_on_plane(&circle, plane);
                let pcurve_sphere = circle_pcurve_on_sphere(&circle, sphere);
                let (pcurve_on_a, pcurve_on_b) = if plane_is_f1 {
                    (Some(pcurve_plane), Some(pcurve_sphere))
                } else {
                    (Some(pcurve_sphere), Some(pcurve_plane))
                };

                let (u_ax_p, v_ax_p) = crate::inttools::edge_face::plane_local_basis(plane);
                // Find arc endpoints clipped to each face's domain.
                // For the planar face, find where the circle intersects the face boundary edges.
                let plane_fi = if plane_is_f1 { f1 } else { f2 };
                let plane_boundary_pts = self.ds.face_boundary_points(plane_fi);
                // Parametric angles of circle points that intersect the plane face boundary.
                let mut clip_angles: Vec<f64> = Vec::new();
                // Use dense angular sampling to find entry/exit points.
                let n_samples = 64usize;
                for i in 0..n_samples {
                    let theta = (i as f64) / (n_samples as f64) * std::f64::consts::TAU;
                    let pt = circle.center + circle.radius * (theta.cos() * u_ax_p + theta.sin() * v_ax_p);
                    let inside = crate::inttools::edge_face::point_in_planar_face_with_tol(
                        pt, plane, &plane_boundary_pts, crate::tolerance::TOLERANCE_ABS * 1000.0);
                    let next_i = (i + 1) % n_samples;
                    let next_theta = (next_i as f64) / (n_samples as f64) * std::f64::consts::TAU;
                    let next_pt = circle.center + circle.radius * (next_theta.cos() * u_ax_p + next_theta.sin() * v_ax_p);
                    let next_inside = crate::inttools::edge_face::point_in_planar_face_with_tol(
                        next_pt, plane, &plane_boundary_pts, crate::tolerance::TOLERANCE_ABS * 1000.0);
                    if inside != next_inside {
                        // Bisect to find precise boundary crossing
                        let t0 = theta;
                        let t1 = if next_theta <= theta { next_theta + std::f64::consts::TAU } else { next_theta };
                        let mut lo = t0;
                        let mut hi = t1;
                        for _ in 0..12 {
                            let mid = (lo + hi) * 0.5;
                            let mp = circle.center + circle.radius * (mid.cos() * u_ax_p + mid.sin() * v_ax_p);
                            let m_inside = crate::inttools::edge_face::point_in_planar_face_with_tol(
                                mp, plane, &plane_boundary_pts, crate::tolerance::TOLERANCE_ABS * 1000.0);
                            if m_inside == inside { lo = mid; } else { hi = mid; }
                        }
                        let cross_angle = (lo + hi) * 0.5;
                        if cross_angle >= 0.0 && cross_angle < std::f64::consts::TAU {
                            clip_angles.push(cross_angle);
                        }
                    }
                }
                clip_angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
                clip_angles.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
                if cfg!(debug_assertions) && std::env::var("RCAD_DEBUG_FF").is_ok() {
                    eprintln!("[DBG_FF] plane-sphere clip: {} angles from {} samples",
                        clip_angles.len(), n_samples);
                }

                let (v_start, v_end, actual_t0, actual_t1) = if clip_angles.len() >= 2 {
                    // Use the first two distinct clip angles as the arc endpoints
                    // (the face interior portion of the circle).
                    let at0 = clip_angles[0];
                    let at1 = clip_angles[1];
                    let p0 = circle.center + circle.radius * (at0.cos() * u_ax_p + at0.sin() * v_ax_p);
                    let p1 = circle.center + circle.radius * (at1.cos() * u_ax_p + at1.sin() * v_ax_p);
                    const IC_VERTEX_MERGE_TOL: f64 = 1e-2;
                    let va = self.ds.find_vertex_near(p0, IC_VERTEX_MERGE_TOL)
                        .unwrap_or_else(|| self.ds.add_vertex(p0));
                    let vb = self.ds.find_vertex_near(p1, IC_VERTEX_MERGE_TOL)
                        .unwrap_or_else(|| self.ds.add_vertex(p1));
                    (va, vb, at0, at1)
                } else {
                    // No valid clip found — fallback to full circle
                    let p_start = circle.center + circle.radius * (effective_t0.cos() * u_ax_p + effective_t0.sin() * v_ax_p);
                    let p_end = circle.center + circle.radius * (effective_t1.cos() * u_ax_p + effective_t1.sin() * v_ax_p);
                    const IC_VERTEX_MERGE_TOL: f64 = 1e-2;
                    let is_closed = p_start.distance_squared(p_end) < TOLERANCE_ABS_SQ;
                    let (v_start, v_end) = if is_closed {
                        let v = self.ds.find_vertex_near(p_start, IC_VERTEX_MERGE_TOL)
                            .unwrap_or_else(|| self.ds.add_vertex(p_start));
                        (v, v)
                    } else {
                        (self.ds.find_vertex_near(p_start, IC_VERTEX_MERGE_TOL)
                            .unwrap_or_else(|| self.ds.add_vertex(p_start)),
                         self.ds.find_vertex_near(p_end, IC_VERTEX_MERGE_TOL)
                            .unwrap_or_else(|| self.ds.add_vertex(p_end)))
                    };
                    (v_start, v_end, effective_t0, effective_t1)
                };
                // OCCT-aligned: inherit tolerance from parent faces (BRep_Tool::Tolerance).
                // OCCT vertices on edges carry source edge/face tolerances (typically 1e-4 to 1e-6);
                // rcad defaults to TOLERANCE_ABS (1e-7) which is too tight for pcurve comparison.
                let parent_tol = self.ds.faces[f1].geom_tol
                    .max(self.ds.faces[f2].geom_tol)
                    .max(self.seam_shift_tol);
                if v_start < self.ds.vertices.len() {
                    self.ds.vertices[v_start].geom_tol = self.ds.vertices[v_start].geom_tol.max(parent_tol);
                }
                if v_end < self.ds.vertices.len() {
                    self.ds.vertices[v_end].geom_tol = self.ds.vertices[v_end].geom_tol.max(parent_tol);
                }

                let curve_idx = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [actual_t0, actual_t1],
                    pcurve_on_a,
                    pcurve_on_b,
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
                if cfg!(debug_assertions) && std::env::var("RCAD_DEBUG_FF").is_ok() {
                    let ic = &self.ds.intersection_curves[curve_idx];
                    eprintln!("[DBG_FF] IC[{}]: f1={} f2={} start_vertex={} end_vertex={} t_range=[{:.6},{:.6}] n_pbs={}",
                        curve_idx, f1, f2, ic.start_vertex, ic.end_vertex, ic.t_range[0], ic.t_range[1], ic.pave_blocks.len());
                }
            }
        }
    }

    // 鈹€鈹€ Sphere �?Sphere analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    /// OCCT: sphere-sphere intersection
    /// OCCT: sphere-sphere intersection
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

    /// OCCT: sphere-cylinder intersection
    /// OCCT: sphere-cylinder intersection
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

    /// OCCT: cylinder-cylinder intersection
    /// OCCT: cylinder-cylinder intersection
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

    /// OCCT: cylinder face V-range
    /// OCCT: cylinder face V-range
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

    /// OCCT: plane-cylinder intersection
    /// OCCT: plane-cylinder intersection
    pub(crate) fn intersect_plane_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, circle_pcurve_on_plane, ellipse_pcurve_on_cylinder,
            ellipse_pcurve_on_plane, line_pcurve_on_cylinder,
            line_pcurve_on_plane,
        };
        use inttools::plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder};
        use rcad_kernel::CurveEval;
        use std::f64::consts::TAU;

        let result = intersect_plane_cylinder(plane, cyl);

        // Determine which face carries the plane (for correct pcurve_on_a/b assignment)
        let plane_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Plane(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if plane_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        let add_curve = |ds: &mut DS,
                         curve: Curve3,
                         t_range: [f64; 2],
                         pcurve_on_a: Option<Curve2d>,
                         pcurve_on_b: Option<Curve2d>,
                         f1: usize,
                         f2: usize|
         -> usize {
            let p_start = curve.point_at(t_range[0]);
            let p_end = curve.point_at(t_range[1]);
            let v_start = ds.add_vertex(p_start);
            let v_end = ds.add_vertex(p_end);
            let curve_idx = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve,
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
            ds.faces[f1].face_info.curves_sc.insert(curve_idx);
            ds.faces[f2].face_info.curves_sc.insert(curve_idx);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            // �?OCCT-aligned:Propagate IC vertices to all faces sharing boundary edges
            //    (BOPDS_FaceInfo::AppendBlock equivalent).
            propagate_ic_vertices_to_shared_faces(ds, &[v_start, v_end], &[f1, f2]);
            curve_idx
        };

        let mut curve_indices = Vec::new();

        match result {
            PlaneCylinderResult::NoIntersection => return,
            PlaneCylinderResult::TangentLine(line) => {
                // �?OCCT aligned: tangent lines are also valid intersection curves,
                //    used to split the cylinder face. OCCT IntTools_FaceFace::MakeCurve
                //    creates BRep edges for tangent lines too.
                // Clip to the cylinder face's parametric V range along the axis
                // so the intersection curve doesn't extend beyond the actual face.
                let cyl_fi = if plane_is_f1 { f2 } else { f1 };
                let v_range = self.cylinder_face_v_range(cyl_fi, cyl);
                let (pca, pcb) = make_pcurves(
                    line_pcurve_on_plane(&line, plane),
                    line_pcurve_on_cylinder(&line, cyl),
                );
                let ci = add_curve(self.ds, Curve3::Line(line), v_range, pca, pcb, f1, f2);
                curve_indices.push(ci);
            }
            PlaneCylinderResult::TwoLines(l1, l2) => {
                // Clip each line to the cylinder face's parametric V range
                let cyl_fi = if plane_is_f1 { f2 } else { f1 };
                let v_range = self.cylinder_face_v_range(cyl_fi, cyl);
                let (pca1, pcb1) = make_pcurves(
                    line_pcurve_on_plane(&l1, plane),
                    line_pcurve_on_cylinder(&l1, cyl),
                );
                let ci1 = add_curve(
                    self.ds,
                    Curve3::Line(l1),
                    v_range,
                    pca1,
                    pcb1,
                    f1,
                    f2,
                );
                let (pca2, pcb2) = make_pcurves(
                    line_pcurve_on_plane(&l2, plane),
                    line_pcurve_on_cylinder(&l2, cyl),
                );
                let ci2 = add_curve(
                    self.ds,
                    Curve3::Line(l2),
                    v_range,
                    pca2,
                    pcb2,
                    f1,
                    f2,
                );
                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }
            PlaneCylinderResult::Circle(circle) => {
                let (pca, pcb) = make_pcurves(
                    circle_pcurve_on_plane(&circle, plane),
                    circle_pcurve_on_cylinder(&circle, cyl),
                );
                let ci = add_curve(
                    self.ds,
                    Curve3::Circle(circle),
                    [0.0, TAU],
                    pca,
                    pcb,
                    f1,
                    f2,
                );
                curve_indices.push(ci);
                // �?OCCT-aligned: when circle lies on cylinder V-boundary, remove IC from cylinder face.
                //    OCCT PerformLoops does not create sub-faces outside face domain, but rcad's BooleanBuilder does.
                //    Plane face retains the IC (split arc segments) for correct box face splitting.
                let cyl_fi = if plane_is_f1 { f2 } else { f1 };
                let v_range = self.cylinder_face_v_range(cyl_fi, cyl);
                let v = (circle.center - cyl.origin).dot(cyl.axis.normalize());
                let boundary_tol = TOLERANCE_ABS * 1000.0;
                if (v - v_range[0]).abs() < boundary_tol
                    || (v - v_range[1]).abs() < boundary_tol
                {
                    self.ds.faces[cyl_fi].face_info.curves_sc.remove(&ci);
                }
            }
            PlaneCylinderResult::Ellipse(ellipse) => {
                let pca_plane = ellipse_pcurve_on_plane(&ellipse, plane);
                let pcb_cyl = ellipse_pcurve_on_cylinder(&ellipse, cyl);
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cyl);
                let ci = add_curve(
                    self.ds,
                    Curve3::Ellipse(ellipse),
                    [0.0, TAU],
                    pca,
                    pcb,
                    f1,
                    f2,
                );
                curve_indices.push(ci);
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

    pub(crate) fn trim_curve_to_faces(
        ds: &DS,
        curve: &Curve3,
        search_range: [f64; 2],
        f1: usize,
        f2: usize,
    ) -> Option<[f64; 2]> {
        use crate::medial_axis::point_in_polygon_2d;
        use rcad_kernel::projection::closest_point_on_surface;
        use std::f64::consts::TAU;

        const N: usize = 256;

        let face1 = &ds.faces[f1];
        let face2 = &ds.faces[f2];
        let uv_bnd1 = face1.uv_boundary.as_ref()?;
        let uv_bnd2 = face2.uv_boundary.as_ref()?;
        let s1 = &face1.surface;
        let s2 = &face2.surface;

        // UV from 3D point on a surface, normalising u �?[0, 2蟺].
        let uv_on_surface = |surface: &Surface3, p: DVec3| -> DVec2 {
            match surface {
                Surface3::Cone(cone) => {
                    let uv = cone.world_to_uv(p);
                    DVec2::new(if uv.x < 0.0 { uv.x + TAU } else { uv.x }, uv.y)
                }
                Surface3::Sphere(sph) => sph.world_to_uv(p),
                Surface3::Cylinder(cyl) => {
                    let x_ax = cyl.ref_dir.normalize();
                    let y_ax = cyl.axis.cross(x_ax).normalize();
                    let local = p - cyl.origin;
                    let u = local.dot(y_ax).atan2(local.dot(x_ax));
                    DVec2::new(if u < 0.0 { u + TAU } else { u }, local.dot(cyl.axis))
                }
                _ => {
                    let proj = closest_point_on_surface(surface, p, 16);
                    DVec2::new(proj.params.0, proj.params.1)
                }
            }
        };

        // True when the curve point at t is inside *both* faces' UV boundaries.
        // For planar faces the 3D point must actually lie on the plane (not just
        // project there), else an off-surface point would be a false positive.
        let point_in_both = |t: f64| -> bool {
            let pt = curve.point_at(t);
            for (sf, bnd) in &[(s1, uv_bnd1), (s2, uv_bnd2)] {
                if let Surface3::Plane(pl) = sf {
                    if (pt - pl.origin).dot(pl.normal).abs() > TOLERANCE_COORD_SUB {
                        return false;
                    }
                }
                let uv = uv_on_surface(sf, pt);
                if !point_in_polygon_2d(uv, bnd) {
                    return false;
                }
            }
            true
        };

        let [t0, t1] = search_range;
        let step = (t1 - t0) / N as f64;
        let mut seg_start: Option<(usize, f64)> = None;
        let mut segments: Vec<(usize, usize, f64, f64)> = Vec::new();

        for i in 0..=N {
            let t = t0 + step * i as f64;
            let inside = point_in_both(t);

            if inside {
                if seg_start.is_none() {
                    seg_start = Some((i, t));
                }
            } else if let Some((si, st)) = seg_start.take() {
                if t - st > TOLERANCE_LINEAR_ULTRA_STRICT {
                    segments.push((si, i, st, t));
                }
            }
        }
        if let Some((si, st)) = seg_start.take() {
            if t1 - st > TOLERANCE_LINEAR_ULTRA_STRICT {
                segments.push((si, N, st, t1));
            }
        }

        // Longest segment
        let (si, ei, rough_start, rough_end) = segments.into_iter().max_by(|a, b| {
            (a.3 - a.2)
                .partial_cmp(&(b.3 - b.2))
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;

        // 鈹€鈹€ binary-search refinement of both endpoints 鈹€鈹€
        // Start: between sample (si-1, outside) and sample (si, inside).
        let refined_start = if si > 0 {
            let t_out = t0 + step * (si - 1) as f64;
            let mut lo = t_out;   // outside
            let mut hi = rough_start; // inside
            for _ in 0..48 {
                let mid = 0.5 * (lo + hi);
                if point_in_both(mid) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            hi
        } else {
            rough_start
        };

        // End: between sample (ei-1, inside) and sample (ei, outside).
        let refined_end = if ei < N {
            let t_in = t0 + step * (ei - 1) as f64;
            let t_out = t0 + step * ei as f64;
            let mut lo = t_in;  // inside
            let mut hi = t_out; // outside
            for _ in 0..48 {
                let mid = 0.5 * (lo + hi);
                if point_in_both(mid) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        } else {
            rough_end
        };

        if refined_end - refined_start > TOLERANCE_LINEAR_ULTRA_STRICT {
            Some([refined_start, refined_end])
        } else {
            None
        }
    }

    // 鈹€鈹€ Plane �?Cone analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    /// OCCT: plane-cone intersection
    /// OCCT: plane-cone intersection
    pub(crate) fn intersect_plane_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        cone: &ConicalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cone, circle_pcurve_on_plane, ellipse_pcurve_on_cone,
            ellipse_pcurve_on_plane, fallback_pcurve_by_projection,
            line_pcurve_on_cone, line_pcurve_on_plane, sampled_pcurve_on_cone,
        };
        use inttools::plane_cone::{PlaneConicalResult, intersect_plane_cone};
        use std::f64::consts::TAU;

        // Determine which face carries the plane
        let plane_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Plane(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if plane_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        // Helper: push a generic curve and register it with both faces.
        let add_curve = |ds: &mut DS,
                         curve: Curve3,
                         t_range: [f64; 2],
                         pcurve_on_a: Option<Curve2d>,
                         pcurve_on_b: Option<Curve2d>,
                         f1: usize,
                         f2: usize|
         -> usize {
            let p_start = curve.point_at(t_range[0]);
            let p_end = curve.point_at(t_range[1]);
            let v_start = ds.add_vertex(p_start);
            let v_end = ds.add_vertex(p_end);
            let ci = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve,
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
            // �?OCCT-aligned:Propagate IC vertices to all faces sharing boundary edges
            //    (BOPDS_FaceInfo::AppendBlock equivalent).
            propagate_ic_vertices_to_shared_faces(ds, &[v_start, v_end], &[f1, f2]);
            ci
        };

        let mut curve_indices = Vec::new();

        match intersect_plane_cone(plane, cone) {
            PlaneConicalResult::NoIntersection | PlaneConicalResult::Point(_) => {
                return;
            }

            PlaneConicalResult::SingleLine(line) => {
                if let Some(trimmed) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Line(line),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca_plane = line_pcurve_on_plane(&line, plane);
                    let pcb_cone = line_pcurve_on_cone(&line, cone);
                    let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                    let ci = add_curve(
                        self.ds,
                        Curve3::Line(line),
                        trimmed,
                        pca,
                        pcb,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci);
                }
            }

            PlaneConicalResult::TwoLines(l1, l2) => {
                if let Some(t1) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Line(l1),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca1 = line_pcurve_on_plane(&l1, plane);
                    let pcb1 = line_pcurve_on_cone(&l1, cone);
                    let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                    let ci1 = add_curve(
                        self.ds,
                        Curve3::Line(l1),
                        t1,
                        pca1,
                        pcb1,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci1);
                }

                if let Some(t2) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Line(l2),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca2 = line_pcurve_on_plane(&l2, plane);
                    let pcb2 = line_pcurve_on_cone(&l2, cone);
                    let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                    let ci2 = add_curve(
                        self.ds,
                        Curve3::Line(l2),
                        t2,
                        pca2,
                        pcb2,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci2);
                }
            }

            PlaneConicalResult::Circle(circle) => {
                let pca_plane = circle_pcurve_on_plane(&circle, plane);
                let pcb_cone = circle_pcurve_on_cone(&circle, cone);
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci = add_curve(self.ds, Curve3::Circle(circle), [0.0, TAU], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::Ellipse(ellipse) => {
                let pca_plane = ellipse_pcurve_on_plane(&ellipse, plane);
                let pcb_cone = ellipse_pcurve_on_cone(&ellipse, cone);
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci = add_curve(self.ds, Curve3::Ellipse(ellipse), [0.0, TAU], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::Parabola(parabola) => {
                if let Some(trimmed) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Parabola(parabola),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca_plane = fallback_pcurve_by_projection(
                        &Curve3::Parabola(parabola),
                        &trimmed,
                        &Surface3::Plane(*plane),
                    );
                    let pcb_cone = sampled_pcurve_on_cone(
                        &Curve3::Parabola(parabola),
                        &trimmed,
                        cone,
                    );
                    let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                    let ci = add_curve(
                        self.ds,
                        Curve3::Parabola(parabola),
                        trimmed,
                        pca,
                        pcb,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci);
                }
            }

            PlaneConicalResult::Hyperbola(hyperbola) => {
                if let Some(trimmed) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Hyperbola(hyperbola),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca_plane = fallback_pcurve_by_projection(
                        &Curve3::Hyperbola(hyperbola),
                        &trimmed,
                        &Surface3::Plane(*plane),
                    );
                    let pcb_cone = sampled_pcurve_on_cone(
                        &Curve3::Hyperbola(hyperbola),
                        &trimmed,
                        cone,
                    );
                    let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                    let ci = add_curve(
                        self.ds,
                        Curve3::Hyperbola(hyperbola),
                        trimmed,
                        pca,
                        pcb,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci);
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

    // 鈹€鈹€ Cylinder �?Cone analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    /// OCCT: cylinder-cone intersection
    /// OCCT: cylinder-cone intersection
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

    /// OCCT: cone-cone intersection
    /// OCCT: cone-cone intersection
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

    /// OCCT: register torus intersection results
    /// OCCT: register torus intersection results
    pub(crate) fn register_torus_intersection(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
        torus_is_f1: bool,
    ) {
        use inttools::intss::{intersect_surfaces_with_density_tol, SurfaceCurve};
        use inttools::pcurve_derive::polyline_pcurve_by_projection;

        let result = intersect_surfaces_with_density_tol(s1, s2, 48, self.ff_tol(f1, f2));
        if result.is_empty() {
            return;
        }

        for sir in &result.curves {
            match &sir.curve_3d {
                SurfaceCurve::Circle(circle) => {
                    // Only split into half-circles for torus脳cylinder intersections where the
                    // full circle spans 100% of cylinder U (triggers has_full_wrap fallback).
                    // For other surface types the full circle is simpler and more robust.
                    // Note: s1 is always Torus by calling convention, s2 is the other surface.
                    if matches!(s2, Surface3::Cylinder(_)) {
                        let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                            if torus_is_f1 { (Some(a.clone()), Some(b.clone())) }
                            else { (Some(b.clone()), Some(a.clone())) }
                        } else { (None, None) };

                        let mut curve_indices = Vec::new();
                        for (t0, t1) in [(0.0, std::f64::consts::PI), (std::f64::consts::PI, std::f64::consts::TAU)] {
                            let pts = sample_circle_arc(circle, t0, t1, 16);
                            if pts.len() < 2 { continue; }
                            let v_start = self.ds.add_vertex(pts[0]);
                            let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                let curve_idx = self.ds.intersection_curves.len();
                eprintln!("[IC] CREATE ci={} f1={} f2={} sv={} ev={}", curve_idx, f1, f2, v_start, v_end);
                            self.ds.intersection_curves.push(IntersectionCurve {
                                curve: Curve3::Circle(*circle),
                                polyline: vec![],
                                start_vertex: v_start,
                                end_vertex: v_end,
                                t_range: [t0, t1],
                                pcurve_on_a: pca.clone(),
                                pcurve_on_b: pcb.clone(),
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

                            curve_indices.push(curve_idx);
                        }

                        if !curve_indices.is_empty() {
                            self.ds.interferences.push(Interference::FaceFace {
                                f1, f2, curves: curve_indices, points: vec![],
                            });
                        }
                    } else {
                        let pts = sample_circle_arc(circle, 0.0, std::f64::consts::TAU, 32);
                        if pts.len() < 2 { continue; }
                        let v_start = self.ds.add_vertex(pts[0]);
                        let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                        let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                            if torus_is_f1 { (Some(a.clone()), Some(b.clone())) }
                            else { (Some(b.clone()), Some(a.clone())) }
                        } else { (None, None) };

                        let curve_idx = self.ds.intersection_curves.len();
                        self.ds.intersection_curves.push(IntersectionCurve {
                            curve: Curve3::Circle(*circle),
                            polyline: vec![],
                            start_vertex: v_start,
                            end_vertex: v_end,
                            t_range: [0.0, std::f64::consts::TAU],
                            pcurve_on_a: pca,
                            pcurve_on_b: pcb,
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
                            f1, f2, curves: vec![curve_idx], points: vec![],
                        });
                    }
                }
                SurfaceCurve::Polyline(pts) => {
                    if pts.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                    let arc_len: f64 = pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
                    let dir = (pts[pts.len() - 1] - pts[0]).normalize_or_zero();

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (
                            polyline_pcurve_by_projection(pts, s1),
                            polyline_pcurve_by_projection(pts, s2),
                        )
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: pts[0],
                            direction: if dir.length_squared() > 0.5 { dir } else { DVec3::X },
                        }),
                        polyline: pts.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
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
                SurfaceCurve::Ellipse(ellipse) => {
                    let pts = sample_circle_arc(
                        &Circle3 {
                            center: ellipse.center,
                            normal: ellipse.normal,
                            radius: ellipse.major_radius,
                        },
                        0.0,
                        std::f64::consts::TAU,
                        32,
                    );
                    if pts.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (None, None)
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Ellipse(*ellipse),
                        polyline: vec![],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, std::f64::consts::TAU],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
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
                SurfaceCurve::Line(line) => {
                    let pts = self.ds.face_boundary_points(f1);
                    let pts2 = self.ds.face_boundary_points(f2);
                    let bbox1_min = pts.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
                    let bbox1_max = pts.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));
                    let bbox2_min = pts2.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
                    let bbox2_max = pts2.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));

                    let lo = bbox1_min.min(bbox2_min);
                    let hi = bbox1_max.max(bbox2_max);
                    let extent = (hi - lo).length() * 0.5 + 1.0;

                    let p_start = line.origin + line.direction * (-extent);
                    let p_end = line.origin + line.direction * extent;

                    let v_start = self.ds.add_vertex(p_start);
                    let v_end = self.ds.add_vertex(p_end);

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (None, None)
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(*line),
                        polyline: vec![p_start, p_end],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [-extent, extent],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
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
                SurfaceCurve::BSplineCurve(b) => {
                    // Sample the BSpline to produce a polyline for face splitting.
                    use rcad_kernel::geom::CurveEval;
                    let n_samples = 33_usize;
                    let mut pts: Vec<DVec3> = Vec::with_capacity(n_samples);
                    for i in 0..n_samples {
                        let t = i as f64 / (n_samples - 1) as f64;
                        pts.push(b.point_at(t));
                    }
                    if pts.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                    let arc_len: f64 = pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
                    let dir = (pts[pts.len() - 1] - pts[0]).normalize_or_zero();

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (
                            polyline_pcurve_by_projection(&pts, s1),
                            polyline_pcurve_by_projection(&pts, s2),
                        )
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::BSpline((**b).clone()),
                        polyline: pts.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
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
                SurfaceCurve::Point(_) | SurfaceCurve::Parabola(_) | SurfaceCurve::Hyperbola(_) => {
                    // Skip degenerate / unsupported curve types for now
                }
            }
        }
    }

    /// OCCT: torus-plane intersection
    /// OCCT: torus-plane intersection
    pub(crate) fn intersect_torus_plane_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        plane: &Plane,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Plane(*plane);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    /// OCCT: torus-sphere intersection
    /// OCCT: torus-sphere intersection
    pub(crate) fn intersect_torus_sphere_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        sphere: &SphericalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Sphere(*sphere);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    /// OCCT: torus-cylinder intersection
    /// OCCT: torus-cylinder intersection
    pub(crate) fn intersect_torus_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        cylinder: &CylindricalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Cylinder(*cylinder);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    /// OCCT: torus-cone intersection
    /// OCCT: torus-cone intersection
    pub(crate) fn intersect_torus_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        cone: &ConicalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Cone(*cone);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    /// OCCT: torus-torus intersection
    /// OCCT: torus-torus intersection
    pub(crate) fn intersect_torus_torus_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus1: &ToroidalSurface,
        torus2: &ToroidalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus1);
        let s2 = Surface3::Torus(*torus2);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    /// OCCT: sphere-cone intersection
    /// OCCT: sphere-cone intersection
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
}
