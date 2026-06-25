use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::geom::CurveEval;
use crate::bopds::ds::{DS, DSEdge, DSRepOnFace, Interference, IntersectionCurve, ShapeOrigin};
use crate::bopds::pave::*;
use crate::bvh::Bvh;
use crate::tolerance::*;
use crate::inttools;
use crate::inttools::context::Context as IntToolsContext;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::pave_filler::helpers::*;
use rcad_kernel::closest_point_on_curve;
use super::propagate_ic_vertices_to_shared_faces;

impl<'a> super::PaveFiller<'a> {

    pub(crate) fn find_face_face_curve_indices(&self, f1: usize, f2: usize) -> Option<Vec<usize>> {
        for inf in &self.ds.interferences {
            if let Interference::FaceFace { f1: a, f2: b, curves, .. } = inf {
                if *a == f1 && *b == f2 {
                    return Some(curves.clone());
                }
            }
        }
        None
    }

    pub(crate) fn sampled_face_boundary_points(&self, face_idx: usize, samples_per_edge: usize) -> Vec<DVec3> {
        let mut pts = Vec::new();
        for &ei in &self.ds.faces[face_idx].boundary_edges {
            if let Some(edge) = self.ds.edges.get(ei) {
                let [t0, t1] = edge.t_range;
                let n = samples_per_edge.max(1);
                for k in 0..=n {
                    let t = t0 + (t1 - t0) * k as f64 / n as f64;
                    let p = edge.curve.point_at(t);
                    if p.is_finite() {
                        pts.push(p);
                    }
                }
            }
        }
        if pts.is_empty() {
            self.ds.face_boundary_points(face_idx)
        } else {
            pts
        }
    }

    pub(crate) fn closest_point_on_boundary_samples(&self, point: DVec3, samples: &[DVec3]) -> DVec3 {
        samples
            .iter()
            .copied()
            .min_by(|a, b| {
                let da = (*a - point).length_squared();
                let db = (*b - point).length_squared();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(point)
    }

    pub(crate) fn snap_polyline_endpoints_to_face_boundaries(
        &self,
        chain: &mut Vec<DVec3>,
        f1: usize,
        f2: usize,
    ) {
        if chain.len() < 2 {
            return;
        }

        let boundary_a = self.sampled_face_boundary_points(f1, 12);
        let boundary_b = self.sampled_face_boundary_points(f2, 12);
        if boundary_a.is_empty() || boundary_b.is_empty() {
            return;
        }

        let snap_start_a = self.closest_point_on_boundary_samples(chain[0], &boundary_a);
        let snap_start_b = self.closest_point_on_boundary_samples(chain[0], &boundary_b);
        let snap_end_a = self.closest_point_on_boundary_samples(chain[chain.len() - 1], &boundary_a);
        let snap_end_b = self.closest_point_on_boundary_samples(chain[chain.len() - 1], &boundary_b);

        let choose_better = |orig: DVec3, p1: DVec3, p2: DVec3| {
            let d1 = (p1 - orig).length_squared();
            let d2 = (p2 - orig).length_squared();
            if d1 <= d2 { p1 } else { p2 }
        };

        let start = choose_better(chain[0], snap_start_a, snap_start_b);
        let end = choose_better(chain[chain.len() - 1], snap_end_a, snap_end_b);

        // Only snap if it is a local correction rather than a gross relocation.
        let local_scale = chain
            .windows(2)
            .map(|w| (w[1] - w[0]).length())
            .filter(|d| d.is_finite() && *d > 0.0)
            .fold(f64::INFINITY, f64::min)
            .min(1.0);
        let snap_tol = (local_scale * 4.0)
            .max(TOLERANCE_RETRY_LADDER_COARSE)
            .max(self.ff_tol(f1, f2));

        if (start - chain[0]).length() <= snap_tol {
            chain[0] = start;
        }
        if (end - chain[chain.len() - 1]).length() <= snap_tol {
            let last = chain.len() - 1;
            chain[last] = end;
        }
    }

    pub(crate) fn check_seam_edge_shift(&self, f1: usize, f2: usize) -> Option<SeamEdgeShift> {
        let s1 = &self.ds.faces[f1].surface;
        let s2 = &self.ds.faces[f2].surface;

        // Skip if both faces are Planes (seam edges only exist on periodic surfaces)
        if matches!(s1, Surface3::Plane(_)) && matches!(s2, Surface3::Plane(_)) {
            return None;
        }

        for &e1 in &self.ds.faces[f1].boundary_edges {
            let is_closed1 = self.is_seam_edge(e1, f1);
            for &e2 in &self.ds.faces[f2].boundary_edges {
                let is_closed2 = self.is_seam_edge(e2, f2);
                if !is_closed1 && !is_closed2 {
                    continue;
                }

                // Look for EE interference between this edge pair
                for inf in &self.ds.interferences {
                    if let Interference::EdgeEdge {
                        e1: ee1,
                        e2: ee2,
                        point,
                        new_vertex,
                        ..
                    } = inf
                    {
                        if !((*ee1 == e1 && *ee2 == e2) || (*ee1 == e2 && *ee2 == e1)) {
                            continue;
                        }

                        // Project the EE vertex point onto both edges' 3D curves
                        // (OCCT: GeomAPI_ProjectPointOnCurve)
                        let curve1 = &self.ds.edges[e1].curve;
                        let curve2 = &self.ds.edges[e2].curve;
                        let proj1 = closest_point_on_curve(curve1, *point, 64);
                        let proj2 = closest_point_on_curve(curve2, *point, 64);

                        let a_p1 = proj1.point;
                        let a_p2 = proj2.point;
                        let shift_dist = a_p1.distance(a_p2);

                        // OCCT-aligned: the seam edge shift is a SMALL tolerance
                        // correction, not a geometric transformation.  Verify both
                        // projections are close to the EE vertex �?if either is
                        // far, the vertex is not near both edges and shifting would
                        // be invalid (e.g. sphere center jumps by 1 unit).
                        let vtx_pt = *point;
                        let d1 = a_p1.distance(vtx_pt);
                        let d2 = a_p2.distance(vtx_pt);
                        // OCCT's shift is a sub-tolerance adjustment.  A projection
                        // error exceeding 1e-4 means the vertex is not on this edge.
                        let sanity_tol = TOLERANCE_ABS * 1000.0;
                        if d1 > sanity_tol || d2 > sanity_tol {
                            continue;
                        }

                        // Check if the shift exceeds vertex tolerance
                        let vtx_tol = self.ds.vertices[*new_vertex].geom_tol;
                        if shift_dist > vtx_tol {
                            // OCCT: shift the face with the closed/seam edge
                            let shift_vector = if is_closed1 {
                                a_p2 - a_p1 // Shift f1: move aP1 toward aP2
                            } else {
                                a_p1 - a_p2 // Shift f2: move aP2 toward aP1
                            };

                            return Some(SeamEdgeShift {
                                shift_vector,
                                shift_value: shift_dist,
                                shifted_face: if is_closed1 { 1 } else { 2 },
                            });
                        }
                    }
                }
            }
        }
        None
    }

    pub(crate) fn reverse_seam_edge_shift(&mut self, f1: usize, f2: usize, shift: &SeamEdgeShift) {
        let inv_vec = if shift.shifted_face == 1 {
            -shift.shift_vector
        } else {
            shift.shift_vector
        };

        // Collect curve indices from the FaceFace interference for this pair
        let mut curve_indices: Vec<usize> = Vec::new();
        for inf in &self.ds.interferences {
            if let Interference::FaceFace {
                f1: a,
                f2: b,
                curves,
                ..
            } = inf
            {
                if (*a == f1 && *b == f2) || (*a == f2 && *b == f1) {
                    curve_indices = curves.clone();
                    break;
                }
            }
        }

        // Reverse shift on each curve
        for &ci in &curve_indices {
            if ci >= self.ds.intersection_curves.len() {
                continue;
            }
            let ic = &mut self.ds.intersection_curves[ci];

            // Translate 3D curve back by inverse shift
            ic.curve = translate_curve3(&ic.curve, inv_vec);

            // Translate polyline points if any
            for p in &mut ic.polyline {
                *p += inv_vec;
            }

            // Translate vertex positions back
            let sv = ic.start_vertex;
            let ev = ic.end_vertex;
            if sv < self.ds.vertices.len() {
                self.ds.vertices[sv].point += inv_vec;
            }
            if ev < self.ds.vertices.len() {
                self.ds.vertices[ev].point += inv_vec;
            }
        }
    }

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

    pub(crate) fn intersect_plane_sphere_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        sphere: &SphericalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_plane, circle_pcurve_on_sphere, fallback_pcurve_by_projection,
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
                // �?OCCT-aligned: clip Circle3 to planar face polygon boundaries
                // map to TWO vertical meridians in sphere UV space. We must create
                // two separate IntersectionCurves �?one per UV branch �?because
                // a single BSpline pcurve cannot span the atan2-wrap discontinuity.
                let is_great = (circle.center - sphere.center).length_squared() < TOLERANCE_ABS_SQ;
                let axis_dot_normal = sphere
                    .axis
                    .normalize()
                    .dot(plane.normal.normalize())
                    .abs();
                let _passes_poles = is_great && axis_dot_normal < TOLERANCE_ABS;

                // �?OCCT-aligned: all plane-sphere ICs go through clip_circle_to_faces unified path.
                //    add_great_circle_curves is disabled �?double-half-arc branches are rcad's own design,
                //    OCCT IntTools_FaceFace clips ICs directly using PutBoundPaveOnCurve.

                // �?OCCT-aligned: clip Circle3 to planar face polygon boundaries
                //    OCCT IntTools_Curve limits range to face boundary at creation time.
                //    rcad: project Circle3 onto plane face 2D polygon,
                //    intersect to get valid parameter range within the face, use its endpoints
                //    as start/end vertices of the curve.
                let clipped_range = self.clip_circle_to_faces(&circle, f1, f2);
                let clipped = match clipped_range {
                    Some(r) => r,
                    None => { return; }
                };
                if clipped[1] - clipped[0] <= TOLERANCE_ABS { return; }
                let (effective_t0, effective_t1) = (clipped[0], clipped[1]);

                if circle.radius <= TOLERANCE_MESH_LEGACY + TOLERANCE_ABS { return; }
                let valid_arc = clipped_range.map(|r| r[1] - r[0]).unwrap_or(0.0);
                if valid_arc <= TOLERANCE_ABS { return; }

                let pcurve_plane = circle_pcurve_on_plane(&circle, plane);
                let pcurve_sphere = fallback_pcurve_by_projection(
                    &Curve3::Circle(circle),
                    &[effective_t0, effective_t1],
                    &Surface3::Sphere(*sphere),
                );
                let (pcurve_on_a, pcurve_on_b) = if plane_is_f1 {
                    (Some(pcurve_plane), Some(pcurve_sphere))
                } else {
                    (Some(pcurve_sphere), Some(pcurve_plane))
                };

                // �?OCCT-aligned: IC endpoints use plane_local_basis (consistent with clip_circle_to_faces)
                //    circle.point_at uses Circle3.normal's any_perpendicular axis,
                //    which may be opposite to plane_local_basis direction, flipping endpoint positions.
                let (u_ax_p, v_ax_p) = crate::inttools::edge_face::plane_local_basis(plane);
                let p_start = circle.center + circle.radius * (effective_t0.cos() * u_ax_p + effective_t0.sin() * v_ax_p);
                let p_end = circle.center + circle.radius * (effective_t1.cos() * u_ax_p + effective_t1.sin() * v_ax_p);
                if p_start.distance_squared(p_end) < TOLERANCE_ABS_SQ { return; }
                // �?OCCT-aligned: try to reuse existing DS vertex (PutPaveOnCurve).
                //    OCCT's IsVertexOnLine detects boundary vertices ON the curve and
                //    places their DS index into the pave block, so the section edge
                //    shares the same TopoDS_Vertex as the boundary edge.  rcad: find
                //    existing vertex within tolerance; only create new if none found.
                //    The tolerance TOLERANCE_ABS*1000 (1e-4) covers intersection noise.
                const IC_VERTEX_MERGE_TOL: f64 = crate::tolerance::TOLERANCE_ABS * 1000.0;
                let v_start = self.ds.find_vertex_near(p_start, IC_VERTEX_MERGE_TOL)
                    .unwrap_or_else(|| self.ds.add_vertex(p_start));
                let v_end = self.ds.find_vertex_near(p_end, IC_VERTEX_MERGE_TOL)
                    .unwrap_or_else(|| self.ds.add_vertex(p_end));
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
                    t_range: [effective_t0, effective_t1],
                    pcurve_on_a,
                    pcurve_on_b,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
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
        }
    }

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

    pub(crate) fn intersect_ff_by_numeric_intss(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
        grid_n: usize,
    ) {
        use inttools::intss::numeric_intss_with_domains;
        use inttools::pcurve_derive::polyline_pcurve_by_projection;

        // Use face-specific UV domains (set up by DS::setup_uv_boundaries)
        // if available.  For cylinders this encodes the actual face height range,
        // ensuring the intersection polyline endpoints fall *inside* the UV
        // boundary rectangle and can be used to split it.
        let dom1 = self.ds.faces[f1]
            .uv_boundary
            .as_ref()
            .and_then(|uv| {
                if uv.len() >= 3 {
                    let u_min = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                    let u_max = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                    let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                    let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
                        return Some([u_min, u_max, v_min, v_max]);
                    }
                }
                None
            });
        let dom2 = self.ds.faces[f2]
            .uv_boundary
            .as_ref()
            .and_then(|uv| {
                if uv.len() >= 3 {
                    let u_min = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                    let u_max = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                    let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                    let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
                        return Some([u_min, u_max, v_min, v_max]);
                    }
                }
                None
            });

        let result = numeric_intss_with_domains(
            s1,
            s2,
            grid_n,
            dom1,
            dom2,
            Some(self.ff_tol(f1, f2)),
        );
        if result.is_empty() {
            return;
        }

        let mut curve_indices = Vec::new();
        for sir in &result.curves {
            let (mut chain, approx_curve) = match &sir.curve_3d {
                crate::inttools::intss::SurfaceCurve::Polyline(pts) => (pts.clone(), None),
                crate::inttools::intss::SurfaceCurve::BSplineCurve(bs) => {
                    // Sample BSpline back to polyline for face-boundary snapping
                    let n = 64usize;
                    let pts: Vec<DVec3> = (0..=n)
                        .map(|i| {
                            let t = i as f64 / n as f64;
                            bs.point_at(t)
                        })
                        .collect();
                    (pts, Some(Curve3::BSpline((**bs).clone())))
                }
                _ => continue,
            };
            if chain.len() < 2 {
                continue;
            }

            self.snap_polyline_endpoints_to_face_boundaries(&mut chain, f1, f2);

            let v_start = self.ds.add_vertex(chain[0]);
            let v_end = self.ds.add_vertex(chain[chain.len() - 1]);

            let arc_len: f64 = chain.windows(2).map(|w| (w[1] - w[0]).length()).sum();
            let dir = (chain[chain.len() - 1] - chain[0]).normalize_or_zero();
            let pcurve_a = sir
                .pcurve_on_a
                .clone()
                .or_else(|| polyline_pcurve_by_projection(&chain, s1));
            let pcurve_b = sir
                .pcurve_on_b
                .clone()
                .or_else(|| polyline_pcurve_by_projection(&chain, s2));

            let curve_idx = self.ds.intersection_curves.len();
            self.ds.intersection_curves.push(IntersectionCurve {
                curve: approx_curve.unwrap_or(Curve3::Line(Line3 {
                    origin: chain[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                })),
                polyline: chain,
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
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
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
            });
        }
    }

    pub(crate) fn intersect_ff_by_marching(&mut self, f1: usize, f2: usize) {
        use inttools::marching::{adaptive_sampling_density, MarchingConfig};

        let s1 = self.ds.faces[f1].surface.clone();
        let s2 = self.ds.faces[f2].surface.clone();

        // OCCT-aligned: use sign-change grid marching (IntTools_FaceFace / IntPatch_ImpPrmIntersection).
        // No BSpline demotion �?BSpline surfaces stay as parametric (ts=0) and use UV grid marching.
        let any_curved = !matches!(&s1, Surface3::Plane(_)) || !matches!(&s2, Surface3::Plane(_));
        if any_curved {
            let char_len = |s: &Surface3| -> f64 {
                match s {
                    Surface3::Sphere(sp) => sp.radius,
                    Surface3::Cylinder(cy) => cy.radius,
                    Surface3::Cone(co) => co.radius.max(0.5),
                    Surface3::Torus(to) => to.major_radius.max(to.minor_radius),
                    Surface3::BSpline(bsp) => {
                        if bsp.control_points.is_empty() {
                            1.0
                        } else {
                            let mut mn = DVec3::splat(f64::INFINITY);
                            let mut mx = DVec3::splat(f64::NEG_INFINITY);
                            for row in &bsp.control_points {
                                for p in row {
                                    mn = mn.min(*p); mx = mx.max(*p);
                                }
                            }
                            (mx - mn).length().max(0.5) * 0.5
                        }
                    }
                    Surface3::Bezier(bez) => {
                        if bez.control_points.is_empty() {
                            1.0
                        } else {
                            let mut mn = DVec3::splat(f64::INFINITY);
                            let mut mx = DVec3::splat(f64::NEG_INFINITY);
                            for row in &bez.control_points {
                                for p in row {
                                    mn = mn.min(*p); mx = mx.max(*p);
                                }
                            }
                            (mx - mn).length().max(0.5) * 0.5
                        }
                    }
                    _ => 1.0,
                }
            };
            let avg_len = (char_len(&s1) + char_len(&s2)) * 0.5;
            let mut grid_n = ((avg_len * 10.0) as usize).max(64).min(256);

            let skew_factor = match (&s1, &s2) {
                (Surface3::Cylinder(c1), Surface3::Cone(c2))
                | (Surface3::Cone(c2), Surface3::Cylinder(c1)) => {
                    let a1 = c1.axis.normalize();
                    let a2 = c2.axis.normalize();
                    let sin_angle = a1.cross(a2).length();
                    (1.0 + sin_angle * 3.0).min(3.0)
                }
                _ => 1.0,
            };
            grid_n = ((grid_n as f64 * skew_factor) as usize).min(256);

            self.intersect_ff_by_numeric_intss(f1, f2, &s1, &s2, grid_n);
            return;
        }

        // Use adaptive sampling density based on surface geometry
        let base_density = 16usize;
        let sampling1 = adaptive_sampling_density(&s1, base_density);
        let sampling2 = adaptive_sampling_density(&s2, base_density);
        // Use the higher density to ensure we don't miss narrow intersections
        let n_u = sampling1.n_u.max(sampling2.n_u);
        let n_v = sampling1.n_v.max(sampling2.n_v);

        let _samples = self.generate_surface_samples_grid(&s1, n_u, n_v);
        // Use multi-scale seed detection for improved robustness
        // Scales: coarse (8x8), medium (16x16), fine (32x32), ultra (64x64)
        let base_step = self.estimate_step_size(&s1, &s2);
        let seed_dedup = (base_step * 2.0).max(self.ff_tol(f1, f2) * 2.0);
        let seeds = inttools::marching::find_seed_points_multiscale(
            &s1,
            &s2,
            |nu, nv| self.generate_surface_samples_grid(&s1, nu, nv),
            &[8, 16, 32, 64],
            seed_dedup,
        );

        if seeds.is_empty() {
            return;
        }

        // Compute a finite bounding box that contains both faces' intersection region.
        // Use boundary vertices (actual face extent) with a generous margin.
        let bounds_from_face = |face_idx: usize| -> (DVec3, DVec3) {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            // Use boundary vertices (from wire edges)
            for &vi in &self.ds.faces[face_idx].boundary_verts {
                let p = self.ds.vertices[vi].point;
                mn = mn.min(p);
                mx = mx.max(p);
            }
            // Also sample boundary edges for curved edges (e.g. circles)
            for &ei in &self.ds.faces[face_idx].boundary_edges {
                if let Some(edge) = self.ds.edges.get(ei) {
                    let [t0, t1] = edge.t_range;
                    for k in 0..=8usize {
                        let t = t0 + (t1 - t0) * k as f64 / 8.0;
                        let p = edge.curve.point_at(t);
                        if p.is_finite() {
                            mn = mn.min(p);
                            mx = mx.max(p);
                        }
                    }
                }
            }
            // If still infinite, use a generous default
            if !mn.is_finite() || !mx.is_finite() {
                mn = DVec3::splat(-10.0);
                mx = DVec3::splat(10.0);
            }
            (mn, mx)
        };

        let (mn1, mx1) = bounds_from_face(f1);
        let (mn2, mx2) = bounds_from_face(f2);
        let margin = 1.0;
        let aabb_min = mn1.min(mn2) - DVec3::splat(margin);
        let aabb_max = mx1.max(mx2) + DVec3::splat(margin);

        // Use adaptive step size based on characteristic lengths
        let char_len = sampling1.characteristic_length.min(sampling2.characteristic_length);
        let step_size = base_step.min(char_len * 0.5).max(TOLERANCE_MESH_LEGACY);

        // Configure marching with convergence monitoring
        let marching_config = MarchingConfig {
            step_size,
            min_step_size: step_size * 0.01,
            max_steps: 500,
            max_oscillations: 3,
            step_reduction_factor: 0.5,
            deflection_tol: step_size * 0.001,
            multiscale_seeds: true,
        };

        let mut curve_indices = Vec::new();
        // Track all points already covered by marched curves, to deduplicate
        // seeds that trace the same intersection curve.
        let mut covered_points: Vec<DVec3> = Vec::new();
        let ff = self.ff_tol(f1, f2);
        let dedup_tol = (step_size * 3.0).max(ff * 2.0);

        for seed in seeds {
            // Skip if this seed is near any point already covered by a previous curve
            if covered_points
                .iter()
                .any(|&cp| (cp - seed).length_squared() < dedup_tol * dedup_tol)
            {
                continue;
            }

            let curve = inttools::marching::march_intersection_with_config(
                &s1,
                &s2,
                seed,
                &marching_config,
                |p: DVec3| p.cmpge(aabb_min).all() && p.cmple(aabb_max).all(),
            );

            if curve.points.len() < 2 {
                continue;
            }

            // Mark all curve points as covered (sample every few for efficiency)
            for (i, &p) in curve.points.iter().enumerate() {
                if i % 5 == 0 {
                    covered_points.push(p);
                }
            }

            let v_start = self.ds.add_vertex(curve.points[0]);
            let v_end = self.ds.add_vertex(curve.points[curve.points.len() - 1]);

            let curve_idx = self.ds.intersection_curves.len();
            // Compute arc-length for t_range
            let arc_len: f64 = curve
                .points
                .windows(2)
                .map(|w| (w[1] - w[0]).length())
                .sum();
            let dir = (curve.points[curve.points.len() - 1] - curve.points[0]).normalize_or_zero();
            let t_range = [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)];

            // �?OCCT-aligned:reApprox �?validate pcurves; retry with loose tolerance
            // if validation fails.
            let (pcurve_a, pcurve_b) = self.make_marching_pcurves_with_reapprox(
                &curve.points, &s1, &s2, f1, f2, &t_range,
            );

            // OCCT-aligned: approximate marching polyline to BSpline (MakeCurve / GeomInt_IntSS::MakeBSpline)
            let approx_curve = if curve.points.len() >= 4 {
                crate::inttools::intss::polyline_to_bspline(&curve.points, TOLERANCE_TOL_SCALE_MICRO)
                    .filter(|c| matches!(c, Curve3::BSpline(_)))
            } else {
                None
            };

            self.ds.intersection_curves.push(IntersectionCurve {
                curve: approx_curve.unwrap_or(Curve3::Line(Line3 {
                    origin: curve.points[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                })),
                polyline: curve.points.clone(),
                start_vertex: v_start,
                end_vertex: v_end,
                t_range,
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
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
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
            });
        }
    }

    pub(crate) fn make_marching_pcurves_with_reapprox(
        &self,
        points: &[DVec3],
        s1: &Surface3,
        s2: &Surface3,
        f1: usize,
        f2: usize,
        t_range: &[f64; 2],
    ) -> (Option<Curve2d>, Option<Curve2d>) {
        let uv_bounds1 = s1.default_domain();
        let uv_bounds2 = s2.default_domain();
        let is_u_periodic1 = matches!(s1, Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_));
        let is_u_periodic2 = matches!(s2, Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_));
        let u_per1 = if is_u_periodic1 { Some(std::f64::consts::TAU) } else { None };
        let u_per2 = if is_u_periodic2 { Some(std::f64::consts::TAU) } else { None };

        // Attempt 1: default tolerance
        let pca = inttools::pcurve_derive::polyline_pcurve_by_projection(points, s1);
        let pcb = inttools::pcurve_derive::polyline_pcurve_by_projection(points, s2);

        let valid_a = pca.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::is_curve_valid_2d(pc)
                && inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds1, u_per1, None)
        });
        let valid_b = pcb.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::is_curve_valid_2d(pc)
                && inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds2, u_per2, None)
        });

        if valid_a && valid_b {
            return (pca, pcb);
        }

        // �?OCCT-aligned:reApprox �?fallback with looser validation.
        // Skip the self-intersection check (is_curve_valid_2d) since polyline
        // pcurves from marching can have V-folds that are geometrically correct.
        let valid_a2 = pca.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds1, u_per1, None)
        });
        let valid_b2 = pcb.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds2, u_per2, None)
        });
        if valid_a2 && valid_b2 {
            return (pca, pcb);
        }

        // Final fallback: return pcurves even if invalid �?the builder handles
        // out-of-face pcurves via its own boundary clipping.
        (pca, pcb)
    }

    pub(crate) fn generate_surface_samples(&self, surface: &Surface3, n1: usize, n2: usize) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                inttools::marching::sample_cylinder(cyl, [-20.0, 20.0], n1, n2)
            }
            Surface3::Sphere(sph) => inttools::marching::sample_sphere(sph, n1, n2),
            Surface3::Torus(torus) => inttools::marching::sample_torus(torus, n1, n2),
            Surface3::Plane(plane) => sample_plane(plane, 20.0, n1),
            Surface3::Cone(cone) => sample_cone(cone, 0.01, 20.0, n1, n2),
            // Generic fallback: sample via surface.default_domain() UV grid.
            // Works for BSpline, Bezier, Offset, Revolution, Trimmed, LinearExtrusion.
            _ => sample_surface_generic(surface, n1, n2),
        }
    }

    pub(crate) fn generate_surface_samples_grid(
        &self,
        surface: &Surface3,
        n_u: usize,
        n_v: usize,
    ) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                // u = azimuth index (0..n_u), v = height index (0..n_v)
                // sample_cylinder returns row = height, col = azimuth,
                // so transpose to row = azimuth, col = height for grid indexing.
                // Rebuild in (n_u azimuth) �?(n_v height) order.
                let height_range = [-20.0_f64, 20.0_f64];
                let u_ax = if cyl.axis.x.abs() < 0.9 {
                    cyl.axis.cross(DVec3::X).normalize()
                } else {
                    cyl.axis.cross(DVec3::Y).normalize()
                };
                let v_ax = cyl.axis.cross(u_ax);
                let mut pts = Vec::with_capacity(n_u * n_v);
                for iu in 0..n_u {
                    let theta =
                        2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
                    for iv in 0..n_v {
                        let h = height_range[0]
                            + (height_range[1] - height_range[0]) * iv as f64
                                / (n_v - 1).max(1) as f64;
                        pts.push(
                            cyl.origin
                                + cyl.axis * h
                                + (u_ax * theta.cos() + v_ax * theta.sin()) * cyl.radius,
                        );
                    }
                }
                pts
            }
            Surface3::Sphere(sph) => {
                let u_ax = if sph.axis.x.abs() < 0.9 {
                    sph.axis.cross(DVec3::X).normalize()
                } else {
                    sph.axis.cross(DVec3::Y).normalize()
                };
                let v_ax = sph.axis.cross(u_ax);
                let mut pts = Vec::with_capacity(n_u * n_v);
                for iu in 0..n_u {
                    let theta =
                        2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
                    for iv in 0..n_v {
                        let phi = std::f64::consts::PI * iv as f64 / (n_v - 1).max(1) as f64;
                        pts.push(
                            sph.center
                                + sph.radius
                                    * (sph.axis * phi.cos()
                                        + (u_ax * theta.cos() + v_ax * theta.sin()) * phi.sin()),
                        );
                    }
                }
                pts
            }
            _ => {
                // Fallback: generic UV-grid sampling for BSpline, Bezier, Offset, etc.
                sample_surface_generic(surface, n_u, n_v)
            }
        }
    }

    pub(crate) fn estimate_step_size(&self, s1: &Surface3, s2: &Surface3) -> f64 {
        // Use a fraction of the smallest characteristic dimension
        let size1 = match s1 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Ellipsoid(e) => e.radius_x.min(e.radius_y).min(e.radius_z),
            Surface3::Pipe(p) => p.radius,
            Surface3::Plane(_)
            | Surface3::Helicoid(_)
            | Surface3::BSpline(_)
            | Surface3::TriBezier(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Ruled(_)
            | Surface3::Coons(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        let size2 = match s2 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Ellipsoid(e) => e.radius_x.min(e.radius_y).min(e.radius_z),
            Surface3::Pipe(p) => p.radius,
            Surface3::Plane(_)
            | Surface3::Helicoid(_)
            | Surface3::BSpline(_)
            | Surface3::TriBezier(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Ruled(_)
            | Surface3::Coons(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        size1.min(size2) * 0.1
    }

}
