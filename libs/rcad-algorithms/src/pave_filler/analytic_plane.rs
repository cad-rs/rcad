use super::*;

impl<'a> super::PaveFiller<'a> {
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
                // OCCT-aligned: IntPatch_ImpImpIntersection returns the full analytic
                // curve without pre-clipping to face boundaries.  OCCT PaveFiller clips
                // via PutBoundPaveOnCurve during MakeBlocks (not at curve creation time).
                if circle.radius <= TOLERANCE_MESH_LEGACY + TOLERANCE_ABS { return; }

                let pcurve_plane = circle_pcurve_on_plane(&circle, plane);
                let pcurve_sphere = circle_pcurve_on_sphere(&circle, sphere);
                let (pcurve_on_a, pcurve_on_b) = if plane_is_f1 {
                    (Some(pcurve_plane), Some(pcurve_sphere))
                } else {
                    (Some(pcurve_sphere), Some(pcurve_plane))
                };

                let (u_ax_p, v_ax_p) = crate::inttools::edge_face::plane_local_basis(plane);
                // Full analytic circle: start_vertex == end_vertex, t_range = [0, 2π].
                // OCCT makes a closed TopoDS_Edge (Edge() == aES); rcad uses a single
                // vertex for both ends.  PrepareLines3D will split the closed curve later.
                let p_start = circle.center + circle.radius * (0.0_f64.cos() * u_ax_p + 0.0_f64.sin() * v_ax_p);
                const IC_VERTEX_MERGE_TOL: f64 = 1e-2;
                let v = self.ds.find_vertex_near(p_start, IC_VERTEX_MERGE_TOL)
                    .unwrap_or_else(|| self.ds.add_vertex(p_start));
                let (v_start, v_end) = (v, v);
                let (actual_t0, actual_t1) = (0.0_f64, std::f64::consts::TAU);

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


}
