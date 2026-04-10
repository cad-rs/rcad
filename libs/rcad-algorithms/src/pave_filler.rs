use glam::DVec3;
use rcad_kernel::geom::*;

use crate::bopds::ds::*;
use crate::bopds::pave::*;
use crate::inttools;
use crate::tolerance::*;

/// PaveFiller executes the six intersection passes (OCCT: BOPAlgo_PaveFiller).
pub struct PaveFiller<'a> {
    pub ds: &'a mut DS,
}

impl<'a> PaveFiller<'a> {
    pub fn new(ds: &'a mut DS) -> Self {
        Self { ds }
    }

    /// Execute all intersection passes.
    pub fn perform(&mut self) {
        self.perform_vv();
        self.perform_ve();
        self.perform_ee();
        self.perform_vf();
        self.perform_ef();
        self.perform_ff();
        self.build_split_edges();
    }

    // ─── Pass 1: Vertex-Vertex ─────────────────────────────────────────

    fn perform_vv(&mut self) {
        let a_verts: Vec<usize> = self
            .ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeA))
            .map(|(i, _)| i)
            .collect();
        let b_verts: Vec<usize> = self
            .ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeB))
            .map(|(i, _)| i)
            .collect();

        for &ai in &a_verts {
            for &bi in &b_verts {
                if points_coincide(self.ds.vertices[ai].point, self.ds.vertices[bi].point) {
                    self.ds.interferences.push(Interference::VertexVertex {
                        v1: ai,
                        v2: bi,
                        merged_vertex: ai,
                    });
                }
            }
        }
    }

    // ─── Pass 2: Vertex-Edge ───────────────────────────────────────────

    fn perform_ve(&mut self) {
        let a_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeA);
        let b_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeB);

        for &vi in &a_verts {
            for &ei in &b_edges {
                self.check_vertex_edge(vi, ei);
            }
        }

        let b_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeB);
        let a_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeA);

        for &vi in &b_verts {
            for &ei in &a_edges {
                self.check_vertex_edge(vi, ei);
            }
        }
    }

    fn check_vertex_edge(&mut self, vi: usize, ei: usize) {
        let point = self.ds.vertices[vi].point;
        let edge = &self.ds.edges[ei];
        match &edge.curve {
            Curve3::Line(line) => {
                if let Some(t) = inttools::vertex_ops::vertex_on_line(point, line, edge.t_range) {
                    self.ds.interferences.push(Interference::VertexEdge {
                        vertex: vi,
                        edge: ei,
                        param: t,
                    });
                    self.ds.edges[ei].paves.push(Pave {
                        vertex_idx: vi,
                        param: t,
                    });
                }
            }
            Curve3::Circle(circle) => {
                // Check if point lies on the circle arc
                let v = point - circle.center;
                let dist = v.length();
                if (dist - circle.radius).abs() < TOLERANCE_ABS {
                    let on_plane = v.dot(circle.normal).abs() < TOLERANCE_ABS;
                    if on_plane {
                        // Compute angular parameter
                        let u = if circle.normal.x.abs() < 0.9 {
                            circle.normal.cross(DVec3::X).normalize()
                        } else {
                            circle.normal.cross(DVec3::Y).normalize()
                        };
                        let w = circle.normal.cross(u);
                        let theta = w.dot(v).atan2(u.dot(v));
                        let t_range = edge.t_range;
                        if theta >= t_range[0] - TOLERANCE_ABS
                            && theta <= t_range[1] + TOLERANCE_ABS
                        {
                            self.ds.interferences.push(Interference::VertexEdge {
                                vertex: vi,
                                edge: ei,
                                param: theta,
                            });
                            self.ds.edges[ei].paves.push(Pave {
                                vertex_idx: vi,
                                param: theta,
                            });
                        }
                    }
                }
            }
            _ => {}
        }
    }

    // ─── Pass 3: Edge-Edge ─────────────────────────────────────────────

    fn perform_ee(&mut self) {
        let a_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeA);
        let b_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeB);

        for &ae in &a_edges {
            for &be in &b_edges {
                self.check_edge_edge(ae, be);
            }
        }
    }

    fn check_edge_edge(&mut self, e1: usize, e2: usize) {
        let edge1 = &self.ds.edges[e1];
        let edge2 = &self.ds.edges[e2];

        if let (Curve3::Line(l1), Curve3::Line(l2)) = (&edge1.curve, &edge2.curve)
            && let Some((t1, t2, point)) = intersect_line_line(l1, edge1.t_range, l2, edge2.t_range)
        {
            let new_v = self.ds.add_vertex(point);
            self.ds.interferences.push(Interference::EdgeEdge {
                e1,
                e2,
                point,
                param1: t1,
                param2: t2,
                new_vertex: new_v,
            });
            self.ds.edges[e1].paves.push(Pave {
                vertex_idx: new_v,
                param: t1,
            });
            self.ds.edges[e2].paves.push(Pave {
                vertex_idx: new_v,
                param: t2,
            });
        }
    }

    // ─── Pass 4: Vertex-Face ───────────────────────────────────────────

    fn perform_vf(&mut self) {
        let a_verts = self.verts_of(ShapeOrigin::ShapeA);
        let b_faces = self.faces_of(ShapeOrigin::ShapeB);

        for &vi in &a_verts {
            for &fi in &b_faces {
                self.check_vertex_face(vi, fi);
            }
        }

        let b_verts = self.verts_of(ShapeOrigin::ShapeB);
        let a_faces = self.faces_of(ShapeOrigin::ShapeA);

        for &vi in &b_verts {
            for &fi in &a_faces {
                self.check_vertex_face(vi, fi);
            }
        }
    }

    fn check_vertex_face(&mut self, vi: usize, fi: usize) {
        let point = self.ds.vertices[vi].point;
        let face = &self.ds.faces[fi];

        if let Surface3::Plane(plane) = &face.surface
            && inttools::vertex_ops::vertex_on_plane(point, plane)
        {
            let face_verts = self.ds.face_boundary_points(fi);
            if inttools::edge_face::point_in_planar_face(point, plane, &face_verts) {
                self.ds.interferences.push(Interference::VertexFace {
                    vertex: vi,
                    face: fi,
                });
                self.ds.faces[fi].face_info.vertices_on.insert(vi);
            }
        } else {
            // For curved surfaces, use closest-point projection to check if
            // the vertex lies on the surface within tolerance.
            let surface = face.surface.clone();
            if !matches!(surface, Surface3::Plane(_)) {
                let proj =
                    rcad_kernel::projection::closest_point_on_surface(&surface, point, 16);
                if proj.distance < TOLERANCE_ABS {
                    self.ds.interferences.push(Interference::VertexFace {
                        vertex: vi,
                        face: fi,
                    });
                    self.ds.faces[fi].face_info.vertices_on.insert(vi);
                }
            }
        }
    }

    // ─── Pass 5: Edge-Face ─────────────────────────────────────────────

    fn perform_ef(&mut self) {
        let a_edges = self.edges_of(ShapeOrigin::ShapeA);
        let b_faces = self.faces_of(ShapeOrigin::ShapeB);

        for &ei in &a_edges {
            for &fi in &b_faces {
                self.intersect_edge_face(ei, fi);
            }
        }

        let b_edges = self.edges_of(ShapeOrigin::ShapeB);
        let a_faces = self.faces_of(ShapeOrigin::ShapeA);

        for &ei in &b_edges {
            for &fi in &a_faces {
                self.intersect_edge_face(ei, fi);
            }
        }
    }

    fn intersect_edge_face(&mut self, edge_idx: usize, face_idx: usize) {
        let edge_curve = self.ds.edges[edge_idx].curve.clone();
        let edge_t_range = self.ds.edges[edge_idx].t_range;
        let face_surface = self.ds.faces[face_idx].surface.clone();

        // Dispatch based on curve type × surface type
        let hits: Vec<(DVec3, f64)> = match (&edge_curve, &face_surface) {
            (Curve3::Line(line), Surface3::Plane(plane)) => {
                inttools::edge_face::intersect_line_plane(line, edge_t_range, plane)
                    .into_iter()
                    .map(|h| (h.point, h.edge_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Cylinder(cyl)) => {
                inttools::curve_surface::intersect_line_cylinder(line, edge_t_range, cyl)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Sphere(sph)) => {
                inttools::curve_surface::intersect_line_sphere(line, edge_t_range, sph)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Cone(cone)) => {
                inttools::curve_surface::intersect_line_cone(line, edge_t_range, cone)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Plane(plane)) => {
                inttools::curve_surface::intersect_circle_plane(circle, edge_t_range, plane)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Cylinder(cyl)) => {
                inttools::curve_surface::intersect_circle_cylinder(circle, edge_t_range, cyl)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Sphere(sph)) => {
                inttools::curve_surface::intersect_circle_sphere(circle, edge_t_range, sph)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Cone(cone)) => {
                inttools::curve_surface::intersect_circle_cone(circle, edge_t_range, cone)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            _ => vec![],
        };

        for (point, edge_param) in hits {
            // Verify hit is within face boundary (for planar faces)
            let in_face = match &face_surface {
                Surface3::Plane(plane) => {
                    let face_verts = self.ds.face_boundary_points(face_idx);
                    inttools::edge_face::point_in_planar_face(point, plane, &face_verts)
                }
                _ => true,
            };

            if !in_face {
                continue;
            }

            // Skip if point is an edge endpoint
            let sv = self.ds.edges[edge_idx].start_vertex;
            let ev = self.ds.edges[edge_idx].end_vertex;
            if points_coincide(point, self.ds.vertices[sv].point)
                || points_coincide(point, self.ds.vertices[ev].point)
            {
                continue;
            }

            let new_v = self.ds.add_vertex(point);
            self.ds.interferences.push(Interference::EdgeFace {
                edge: edge_idx,
                face: face_idx,
                point,
                edge_param,
                new_vertex: new_v,
            });
            self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
            self.ds.edges[edge_idx].paves.push(Pave {
                vertex_idx: new_v,
                param: edge_param,
            });
        }
    }

    // ─── Pass 6: Face-Face ─────────────────────────────────────────────

    fn perform_ff(&mut self) {
        let a_faces = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces = self.faces_of(ShapeOrigin::ShapeB);

        for &af in &a_faces {
            for &bf in &b_faces {
                self.intersect_face_face(af, bf);
            }
        }
    }

    fn intersect_face_face(&mut self, f1: usize, f2: usize) {
        let s1 = self.ds.faces[f1].surface.clone();
        let s2 = self.ds.faces[f2].surface.clone();

        match (&s1, &s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                self.intersect_plane_plane_faces(f1, f2, p1, p2);
            }
            (Surface3::Plane(pl), Surface3::Sphere(sph))
            | (Surface3::Sphere(sph), Surface3::Plane(pl)) => {
                self.intersect_plane_sphere_faces(f1, f2, pl, sph);
            }
            (Surface3::Plane(pl), Surface3::Cylinder(cyl))
            | (Surface3::Cylinder(cyl), Surface3::Plane(pl)) => {
                self.intersect_plane_cylinder_faces(f1, f2, pl, cyl);
            }
            (Surface3::Sphere(sph1), Surface3::Sphere(sph2)) => {
                let (sph1, sph2) = (*sph1, *sph2);
                self.intersect_sphere_sphere_faces(f1, f2, &sph1, &sph2);
            }
            (Surface3::Sphere(sph), Surface3::Cylinder(cyl))
            | (Surface3::Cylinder(cyl), Surface3::Sphere(sph)) => {
                let (sph, cyl) = (*sph, *cyl);
                self.intersect_sphere_cylinder_faces(f1, f2, &sph, &cyl);
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                let (c1, c2) = (*c1, *c2);
                self.intersect_cylinder_cylinder_faces(f1, f2, &c1, &c2);
            }
            _ => {
                // General case: numerical marching
                self.intersect_ff_by_marching(f1, f2);
            }
        }
    }

    fn intersect_plane_plane_faces(&mut self, f1: usize, f2: usize, p1: &Plane, p2: &Plane) {
        use inttools::pcurve_derive::line_pcurve_on_plane;

        match inttools::plane_plane::intersect_plane_plane(p1, p2) {
            inttools::plane_plane::PlanePlaneResult::Parallel => {}
            inttools::plane_plane::PlanePlaneResult::Coincident => {
                // Coplanar — handled via coplanar analysis
                self.handle_coplanar_faces(f1, f2, p1);
            }
            inttools::plane_plane::PlanePlaneResult::Line(line) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);

                let range1 = inttools::edge_face::clip_line_to_convex_polygon(&line, p1, &verts1);
                let range2 = inttools::edge_face::clip_line_to_convex_polygon(&line, p2, &verts2);

                if let (Some((t1_min, t1_max)), Some((t2_min, t2_max))) = (range1, range2) {
                    let t_min = t1_min.max(t2_min);
                    let t_max = t1_max.min(t2_max);
                    if t_max - t_min < TOLERANCE_ABS {
                        return;
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
                    });

                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![curve_idx],
                        points: vec![],
                    });

                    self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
            }
        }
    }

    fn handle_coplanar_faces(&mut self, f1: usize, f2: usize, plane: &Plane) {
        let verts1 = self.ds.face_boundary_points(f1);
        let verts2 = self.ds.face_boundary_points(f2);

        let result = inttools::coplanar::analyze_coplanar_faces(&verts1, &verts2, plane);

        if !result.overlap.is_empty() {
            // Record as a FaceFace interference with no curves (coplanar overlap)
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: vec![],
                points: vec![],
            });
        }
    }

    // ── Plane × Sphere analytic face-face intersection ─────────────────────────

    fn intersect_plane_sphere_faces(
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

        match intersect_plane_sphere(plane, sphere) {
            PlaneSphereResult::NoIntersection => {}
            PlaneSphereResult::TangentPoint(pt) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);
                if inttools::edge_face::point_in_planar_face(pt, plane, &verts1)
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
                // Sample the circle and clip to both face boundaries
                let pts = sample_circle_arc(&circle, 0.0, std::f64::consts::TAU, 32);
                if pts.len() < 2 {
                    return;
                }

                let pcurve_plane = circle_pcurve_on_plane(&circle, plane);
                // `circle_pcurve_on_sphere` is only analytically correct when the
                // sphere axis is parallel to the cutting plane normal (i.e. the
                // intersection circle is a latitude line in the sphere's UV domain).
                // When the axis is not aligned with the plane normal, we fall back to
                // projection-based sampling — exactly as `intersect_sphere_sphere_faces`
                // already does — to obtain the correct parameter-space curve.
                let axis_dot_normal = sphere
                    .axis
                    .normalize()
                    .dot(plane.normal.normalize())
                    .abs();
                let pcurve_sphere = if (axis_dot_normal - 1.0).abs() < 1e-6 {
                    // Axis is parallel to plane normal → latitude line is exact.
                    circle_pcurve_on_sphere(&circle, sphere)
                } else {
                    // Axis is NOT aligned → use projection fallback.
                    fallback_pcurve_by_projection(
                        &Curve3::Circle(circle),
                        &[0.0, std::f64::consts::TAU],
                        &Surface3::Sphere(*sphere),
                    )
                };
                let (pcurve_on_a, pcurve_on_b) = if plane_is_f1 {
                    (Some(pcurve_plane), Some(pcurve_sphere))
                } else {
                    (Some(pcurve_sphere), Some(pcurve_plane))
                };

                let v_start = self.ds.add_vertex(pts[0]);
                let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                let curve_idx = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, std::f64::consts::TAU],
                    pcurve_on_a,
                    pcurve_on_b,
                });

                self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    // ── Sphere × Sphere analytic face-face intersection ───────────────────────

    fn intersect_sphere_sphere_faces(
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
        if d < 1e-14 || d >= sph1.radius + sph2.radius || d <= (sph1.radius - sph2.radius).abs() {
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
        // Use projection-based PCurves since the circle may not be a latitude line
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
        });

        self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
        self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    // ── Sphere × Cylinder analytic face-face intersection ─────────────────────

    fn intersect_sphere_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sphere: &SphericalSurface,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, circle_pcurve_on_sphere, fallback_pcurve_by_projection,
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
                });
                ds.faces[f1].face_info.curves_in.insert(curve_idx);
                ds.faces[f2].face_info.curves_in.insert(curve_idx);
                ds.faces[f1].face_info.vertices_in.insert(v_start);
                ds.faces[f1].face_info.vertices_in.insert(v_end);
                ds.faces[f2].face_info.vertices_in.insert(v_start);
                ds.faces[f2].face_info.vertices_in.insert(v_end);
                curve_idx
            };

        // Closure to compute pcurves for one intersection circle.
        // Uses `circle_pcurve_on_sphere` when the sphere axis is parallel to
        // the cylinder axis (so the circle is a true latitude line), otherwise
        // falls back to projection sampling.
        let make_circle_pcurves = |circle: &Circle3| -> (Option<Curve2d>, Option<Curve2d>) {
            let axis_dot = sphere.axis.normalize().dot(cyl.axis.normalize()).abs();
            let pcurve_sph = if (axis_dot - 1.0).abs() < 1e-6 {
                circle_pcurve_on_sphere(circle, sphere)
            } else {
                fallback_pcurve_by_projection(
                    &Curve3::Circle(*circle),
                    &[0.0, TAU],
                    &Surface3::Sphere(*sphere),
                )
            };
            let pcurve_cyl = circle_pcurve_on_cylinder(circle, cyl);
            make_pcurves(pcurve_sph, pcurve_cyl)
        };

        match intersect_sphere_cylinder(sphere, cyl) {
            SphereCylinderResult::NoIntersection => return,
            SphereCylinderResult::General => {
                // Fall back to numeric marching for the quartic case.
                self.intersect_ff_by_marching(f1, f2);
                return;
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
        }
    }

    // ── Cylinder × Cylinder analytic face-face intersection ──────────────────

    fn intersect_cylinder_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cyl1: &CylindricalSurface,
        cyl2: &CylindricalSurface,
    ) {
        use inttools::cylinder_cylinder::{CylinderCylinderResult, intersect_cylinder_cylinder};
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, ellipse_pcurve_on_plane, fallback_pcurve_by_projection,
            line_pcurve_on_cylinder, line_pcurve_on_plane,
        };
        use std::f64::consts::TAU;

        // Determine which face is cyl1 (for pcurve_on_a/b ordering)
        let cyl1_is_f1 = {
            if let Surface3::Cylinder(c) = &self.ds.faces[f1].surface {
                (c.origin - cyl1.origin).length_squared() < 1e-10
                    && (c.axis - cyl1.axis).length_squared() < 1e-10
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
                });
                ds.faces[f1].face_info.curves_in.insert(ci);
                ds.faces[f2].face_info.curves_in.insert(ci);
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
            });
            ds.faces[f1].face_info.curves_in.insert(ci);
            ds.faces[f2].face_info.curves_in.insert(ci);
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
            });
            ds.faces[f1].face_info.curves_in.insert(ci);
            ds.faces[f2].face_info.curves_in.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            ci
        };

        let extent = 20.0_f64;
        let mut curve_indices = Vec::new();

        match intersect_cylinder_cylinder(cyl1, cyl2) {
            CylinderCylinderResult::NoIntersection | CylinderCylinderResult::Coaxial => return,

            CylinderCylinderResult::General => {
                // Fall back to numeric marching for skew/oblique cases.
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
                // PCurves for the cylinder surfaces use projection fallback since
                // these circles are not latitude or generator lines.
                let pca1 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb1 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_circle(self.ds, &c1, pca1, pcb1, f1, f2);

                let pca2 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb2 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_circle(self.ds, &c2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::TwoEllipses(e1, e2) => {
                // Perpendicular Steinmetz unequal-radii.
                // Each ellipse lies in a plane; use ellipse_pcurve_on_plane for the
                // plane PCurve and fallback for the cylinder PCurve.
                let plane1 = rcad_kernel::geom::Plane {
                    origin: e1.center,
                    normal: e1.normal,
                };
                let plane2 = rcad_kernel::geom::Plane {
                    origin: e2.center,
                    normal: e2.normal,
                };

                let pca1 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb1 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_ellipse(self.ds, &e1, pca1, pcb1, f1, f2);

                let pca2 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb2 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_ellipse(self.ds, &e2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);

                let _ = (plane1, plane2); // suppress warnings
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

    // ── Plane × Cylinder analytic face-face intersection ──────────────────────

    fn intersect_plane_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, circle_pcurve_on_plane, ellipse_pcurve_on_plane,
            fallback_pcurve_by_projection, line_pcurve_on_cylinder, line_pcurve_on_plane,
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
            });
            ds.faces[f1].face_info.curves_in.insert(curve_idx);
            ds.faces[f2].face_info.curves_in.insert(curve_idx);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            curve_idx
        };

        let mut curve_indices = Vec::new();

        match result {
            PlaneCylinderResult::NoIntersection => return,
            PlaneCylinderResult::TangentLine(_) => return, // zero-area intersection
            PlaneCylinderResult::TwoLines(l1, l2) => {
                // Clip each line to the face bounding-box extent
                let extent = 20.0_f64;
                let (pca1, pcb1) = make_pcurves(
                    line_pcurve_on_plane(&l1, plane),
                    line_pcurve_on_cylinder(&l1, cyl),
                );
                let ci1 = add_curve(
                    self.ds,
                    Curve3::Line(l1),
                    [-extent, extent],
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
                    [-extent, extent],
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
            }
            PlaneCylinderResult::Ellipse(ellipse) => {
                let pca_plane = ellipse_pcurve_on_plane(&ellipse, plane);
                let pcb_cyl = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(ellipse),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl),
                );
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

    /// For curved×curved face pairs, use numeric_intss_with_density (sign-change
    /// edge marching) which returns ordered polylines without the closure/drift
    /// issues of the gradient marcher.
    fn intersect_ff_by_numeric_intss(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
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

        let result = numeric_intss_with_domains(s1, s2, 64, dom1, dom2);
        if result.is_empty() {
            return;
        }

        let mut curve_indices = Vec::new();
        for sir in &result.curves {
            let chain = match &sir.curve_3d {
                crate::inttools::intss::SurfaceCurve::Polyline(pts) => pts.clone(),
                _ => continue,
            };
            if chain.len() < 2 {
                continue;
            }

            let v_start = self.ds.add_vertex(chain[0]);
            let v_end = self.ds.add_vertex(chain[chain.len() - 1]);

            let arc_len: f64 = chain.windows(2).map(|w| (w[1] - w[0]).length()).sum();
            let dir = (chain[chain.len() - 1] - chain[0]).normalize_or_zero();
            let pcurve_a = polyline_pcurve_by_projection(&chain, s1);
            let pcurve_b = polyline_pcurve_by_projection(&chain, s2);

            let curve_idx = self.ds.intersection_curves.len();
            self.ds.intersection_curves.push(IntersectionCurve {
                curve: Curve3::Line(Line3 {
                    origin: chain[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                }),
                polyline: chain,
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, arc_len.max(1e-10)],
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
            });

            self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
            self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    fn intersect_ff_by_marching(&mut self, f1: usize, f2: usize) {
        use inttools::pcurve_derive::polyline_pcurve_by_projection;

        let s1 = self.ds.faces[f1].surface.clone();
        let s2 = self.ds.faces[f2].surface.clone();

        // For curved-curved surface pairs (e.g. Cylinder × Cylinder), use the
        // sign-change grid marching from numeric_intss_with_density, which produces
        // ordered polylines directly without the closure/drift issues of the gradient
        // marcher.
        let both_curved = !matches!(&s1, Surface3::Plane(_)) && !matches!(&s2, Surface3::Plane(_));
        if both_curved {
            self.intersect_ff_by_numeric_intss(f1, f2, &s1, &s2);
            return;
        }

        // For plane-curved pairs: use gradient marching with grid-based seed finder.
        let (n_u, n_v) = (16usize, 16usize);
        let samples = self.generate_surface_samples_grid(&s1, n_u, n_v);
        let seeds = inttools::marching::find_seed_points_grid(&s1, &s2, &samples, n_u, n_v);

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

        // March each seed
        let step_size = self.estimate_step_size(&s1, &s2);
        let mut curve_indices = Vec::new();
        // Track all points already covered by marched curves, to deduplicate
        // seeds that trace the same intersection curve.
        let mut covered_points: Vec<DVec3> = Vec::new();
        let dedup_tol = step_size * 3.0;

        for seed in seeds {
            // Skip if this seed is near any point already covered by a previous curve
            if covered_points
                .iter()
                .any(|&cp| (cp - seed).length_squared() < dedup_tol * dedup_tol)
            {
                continue;
            }

            let curve = inttools::marching::march_intersection(
                &s1,
                &s2,
                seed,
                step_size,
                500,
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
            let pcurve_a = polyline_pcurve_by_projection(&curve.points, &s1);
            let pcurve_b = polyline_pcurve_by_projection(&curve.points, &s2);
            self.ds.intersection_curves.push(IntersectionCurve {
                curve: Curve3::Line(Line3 {
                    origin: curve.points[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                }),
                polyline: curve.points.clone(),
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, arc_len.max(1e-10)],
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
            });

            self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
            self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    fn generate_surface_samples(&self, surface: &Surface3, n1: usize, n2: usize) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                inttools::marching::sample_cylinder(cyl, [-5.0, 5.0], n1, n2)
            }
            Surface3::Sphere(sph) => inttools::marching::sample_sphere(sph, n1, n2),
            Surface3::Torus(torus) => inttools::marching::sample_torus(torus, n1, n2),
            Surface3::Plane(plane) => sample_plane(plane, 5.0, n1),
            Surface3::Cone(cone) => sample_cone(cone, 0.01, 5.0, n1, n2),
            _ => vec![],
        }
    }

    /// Like `generate_surface_samples` but returns a structured `n_u × n_v` grid
    /// (row-major) so callers can use grid-aware adjacency for seed detection.
    fn generate_surface_samples_grid(
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
                // Rebuild in (n_u azimuth) × (n_v height) order.
                let height_range = [-5.0_f64, 5.0_f64];
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
                // Fallback: flat list
                self.generate_surface_samples(surface, n_u, n_v)
            }
        }
    }

    fn estimate_step_size(&self, s1: &Surface3, s2: &Surface3) -> f64 {
        // Use a fraction of the smallest characteristic dimension
        let size1 = match s1 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Plane(_)
            | Surface3::BSpline(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        let size2 = match s2 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Plane(_)
            | Surface3::BSpline(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        size1.min(size2) * 0.1
    }

    // ─── Edge splitting ────────────────────────────────────────────────

    fn build_split_edges(&mut self) {
        for ei in 0..self.ds.edges.len() {
            let edge = &self.ds.edges[ei];
            if edge.paves.is_empty() {
                // No splits — single pave block spanning entire edge
                let pb = PaveBlock::new(
                    ei,
                    Pave {
                        vertex_idx: edge.start_vertex,
                        param: edge.t_range[0],
                    },
                    Pave {
                        vertex_idx: edge.end_vertex,
                        param: edge.t_range[1],
                    },
                );
                self.ds.edges[ei].pave_blocks = vec![pb];
                continue;
            }

            // Collect all paves including endpoints, sort by parameter
            let mut all_paves = vec![
                Pave {
                    vertex_idx: edge.start_vertex,
                    param: edge.t_range[0],
                },
                Pave {
                    vertex_idx: edge.end_vertex,
                    param: edge.t_range[1],
                },
            ];
            all_paves.extend_from_slice(&edge.paves);
            all_paves.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));

            // Deduplicate paves at the same parameter
            all_paves.dedup_by(|a, b| params_equal(a.param, b.param));

            // Create pave blocks between consecutive paves
            let mut blocks = Vec::new();
            for w in all_paves.windows(2) {
                blocks.push(PaveBlock::new(ei, w[0], w[1]));
            }
            self.ds.edges[ei].pave_blocks = blocks;
        }
    }

    // ─── Helpers ───────────────────────────────────────────────────────

    fn verts_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(origin))
            .map(|(i, _)| i)
            .collect()
    }

    fn edges_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }

    fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }
}

/// Intersect two bounded line segments in 3D. Returns (t1, t2, point) if they
/// cross within tolerance.
fn intersect_line_line(
    l1: &Line3,
    r1: [f64; 2],
    l2: &Line3,
    r2: [f64; 2],
) -> Option<(f64, f64, DVec3)> {
    let d1 = l1.direction;
    let d2 = l2.direction;
    let w0 = l1.origin - l2.origin;

    let a = d1.dot(d1);
    let b = d1.dot(d2);
    let c = d2.dot(d2);
    let d = d1.dot(w0);
    let e = d2.dot(w0);

    let denom = a * c - b * b;
    if denom.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        return None; // parallel
    }

    let t1 = (b * e - c * d) / denom;
    let t2 = (a * e - b * d) / denom;

    // Check within ranges
    if t1 < r1[0] - TOLERANCE_ABS
        || t1 > r1[1] + TOLERANCE_ABS
        || t2 < r2[0] - TOLERANCE_ABS
        || t2 > r2[1] + TOLERANCE_ABS
    {
        return None;
    }

    let p1 = l1.origin + d1 * t1;
    let p2 = l2.origin + d2 * t2;

    if !points_coincide(p1, p2) {
        return None; // skew, don't actually intersect
    }

    Some((t1, t2, (p1 + p2) * 0.5))
}

// ── Sampling helpers for marching seed-point generation ──────────────────────

/// Sample a flat plane (infinite) over a 2D square of side `half_extent*2`
/// centred at `plane.origin`.
fn sample_plane(plane: &Plane, half_extent: f64, n: usize) -> Vec<DVec3> {
    let u = rcad_kernel::any_perpendicular(plane.normal);
    let v = plane.normal.cross(u);
    let mut pts = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let su = -half_extent + 2.0 * half_extent * i as f64 / (n - 1).max(1) as f64;
            let sv = -half_extent + 2.0 * half_extent * j as f64 / (n - 1).max(1) as f64;
            pts.push(plane.origin + u * su + v * sv);
        }
    }
    pts
}

/// Sample a cone surface between heights `h_min` and `h_max` along its axis.
fn sample_cone(
    cone: &ConicalSurface,
    h_min: f64,
    h_max: f64,
    n_theta: usize,
    n_h: usize,
) -> Vec<DVec3> {
    let u = rcad_kernel::any_perpendicular(cone.axis);
    let v = cone.axis.cross(u);
    let tan_h = cone.half_angle_rad.tan();
    let mut pts = Vec::with_capacity(n_theta * n_h);
    for ih in 0..n_h {
        let h = h_min + (h_max - h_min) * ih as f64 / (n_h - 1).max(1) as f64;
        let r = h * tan_h;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = cone.apex + cone.axis * h + (u * theta.cos() + v * theta.sin()) * r;
            pts.push(p);
        }
    }
    pts
}

/// Sample `n` points on a circular arc from `t_start` to `t_end`.
fn sample_circle_arc(circle: &Circle3, t_start: f64, t_end: f64, n: usize) -> Vec<DVec3> {
    use rcad_kernel::CurveEval;
    use rcad_kernel::geom::Curve3;
    let curve = Curve3::Circle(*circle);
    (0..n)
        .map(|i| {
            let t = t_start + (t_end - t_start) * i as f64 / (n - 1).max(1) as f64;
            curve.point_at(t)
        })
        .collect()
}

/// Check whether a point lies within the boundary of a sphere-face, defined by
/// the sphere face boundary vertices (used for tangent-point containment check).
fn point_in_sphere_face(pt: DVec3, boundary_verts: &[DVec3], _ds: &DS) -> bool {
    // Simple bounding-box check: the point should be within the convex hull
    // of the boundary vertices on the sphere surface (rough approximation).
    if boundary_verts.is_empty() {
        return true;
    }
    let cx = boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::max) + 1e-9);
    let cy = boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::max) + 1e-9);
    let cz = boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::max) + 1e-9);
    cx.contains(&pt.x) && cy.contains(&pt.y) && cz.contains(&pt.z)
}
